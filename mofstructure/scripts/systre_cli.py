#!/usr/bin/env python3
from __future__ import annotations
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"

import argparse
import csv
import gc
import hashlib
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

from mofstructure.structure import MOFstructure
from mofstructure.filetyper import convert_numpy_types, load_dict_msgpack, save_dict_msgpack


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compute topology using MOFstructure.get_topology() for files or folders."
    )
    parser.add_argument(
        "inputs",
        nargs="*",
        help="One or more files and/or folders. Not required with "
             "--finalise-only, which only merges existing batches.",
    )
    parser.add_argument(
        "--method",
        choices=["sbus", "ligand_cluster", "all_node", "single_node", "all"],
        default="all_node",
        help="Topology method passed to MOFstructure.get_topology(). "
             "Use 'all' to run every method and record them "
             "all in one entry per structure.",
    )
    parser.add_argument(
        "--decimals",
        type=int,
        default=8,
        help="Number of decimal places used in topology hashing.",
    )
    parser.add_argument(
        "--include-edge-centers",
        action="store_true",
        help="Include edge-center comments in the CRYSTAL2 text.",
    )
    parser.add_argument(
        "--fallback-to-input-cgd",
        action="store_true",
        help="Fallback to input CGD text when needed.",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        default=True,
        help="Recurse into folders (default: True).",
    )
    parser.add_argument(
        "--no-recursive",
        action="store_false",
        dest="recursive",
        help="Do not recurse into folders.",
    )
    parser.add_argument(
        "--flush-every",
        type=int,
        default=1000,
        help="Flush one MessagePack batch every N processed files (default: 100).",
    )
    parser.add_argument(
        "--batch-dir",
        default="Topology_msgpack_batches",
        help="Directory for intermediate MessagePack batch files.",
    )
    parser.add_argument(
        "--msgpack",
        dest="msgpack_out",
        default="Topology.msgpack",
        help="Final merged MessagePack file.",
    )
    parser.add_argument(
        "--json",
        dest="json_out",
        default="Topology.json",
        help="Final merged JSON file.",
    )
    parser.add_argument(
        "--csv",
        dest="csv_out",
        default="Topology.csv",
        help="Final merged CSV summary file.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        default=True,
        help="Resume from existing batch files (default: True).",
    )
    parser.add_argument(
        "--no-resume",
        action="store_false",
        dest="resume",
        help="Ignore existing batch files and recompute everything.",
    )
    parser.add_argument(
        "--finalise-only",
        action="store_true",
        help="Do not compute new files; only merge existing batch files into final outputs.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print full topology dictionaries during processing.",
    )
    return parser


def _collect_input_files(inputs: Sequence[str], recursive: bool = True) -> List[Path]:
    """
    Fast collection: do not pre-read with ASE. Just gather files and let compute fail naturally.
    """
    collected: List[Path] = []
    seen: Set[str] = set()

    for item in inputs:
        path = Path(item)

        if not path.exists():
            print(f"Warning: input not found: {path}")
            continue

        if path.is_file():
            resolved = str(path.resolve())
            if resolved not in seen:
                collected.append(path)
                seen.add(resolved)
            continue

        if path.is_dir():
            iterator = path.rglob("*") if recursive else path.glob("*")
            for subpath in iterator:
                if not subpath.is_file():
                    continue
                resolved = str(subpath.resolve())
                if resolved not in seen:
                    collected.append(subpath)
                    seen.add(resolved)

    return collected


def _safe_stem(path: Path) -> str:
    """
    Create a filesystem- and JSON-friendly name stem.
    """
    stem = path.stem.strip()
    stem = re.sub(r"\s+", "_", stem)
    stem = re.sub(r"[^A-Za-z0-9_.-]", "_", stem)
    return stem or "structure"


def _make_record_key(path: Path) -> str:
    """
    Deterministic unique key based on resolved path.
    This avoids collisions and makes resume safe.
    """
    resolved = str(path.resolve())
    digest = hashlib.sha1(resolved.encode("utf-8")).hexdigest()[:12]
    return f"{_safe_stem(path)}__{digest}"


# the individual node definitions run when method="all"
ALL_METHODS = ("sbus", "ligand_cluster", "all_node", "single_node")

# stand-in record for a node definition that could not be computed
_ERROR_TOPOLOGY = {
    "topology": "ERROR",
    "dimension": None,
    "td10": None,
    "topology_hash": None,
    "cgd": None,
}


def _compute_topology(
    path: Path,
    *,
    method: str,
    decimals: int,
    include_edge_centers: bool,
    fallback_to_input_cgd: bool,
) -> Dict[str, object]:
    '''
    Compute the topology of one structure.

    For a single method the get_topology fields are returned flat. For
    method="all" each node definition in ALL_METHODS is computed and the
    results are nested under a "topologies" mapping keyed by method name, so one
    record carries every net and drops straight into a database.
    '''
    structure = MOFstructure(filename=str(path))

    if method != "all":
        data = structure.get_topology(
            method=method,
            decimals=decimals,
            include_edge_centers=include_edge_centers,
            fallback_to_input_cgd=fallback_to_input_cgd,
        )
        return convert_numpy_types(data)

    topologies = {}
    for node_method in ALL_METHODS:
        try:
            data = structure.get_topology(
                method=node_method,
                decimals=decimals,
                include_edge_centers=include_edge_centers,
                fallback_to_input_cgd=fallback_to_input_cgd,
            )
            topologies[node_method] = convert_numpy_types(data)
        except Exception:
            topologies[node_method] = dict(_ERROR_TOPOLOGY)
    return {"topologies": topologies}


def _write_batch(
    batch_dir: Path,
    batch_index: int,
    batch_data: Dict[str, Dict[str, object]],
    batch_sources: List[str],
) -> Tuple[Path, Path]:
    """
    Write one MessagePack batch plus a tiny .sources checkpoint file.
    """
    batch_dir.mkdir(parents=True, exist_ok=True)

    msgpack_path = batch_dir / f"batch_{batch_index:06d}.msgpack"
    sources_path = batch_dir / f"batch_{batch_index:06d}.sources"

    save_dict_msgpack(batch_data, str(msgpack_path))

    with sources_path.open("w", encoding="utf-8") as handle:
        for source in batch_sources:
            handle.write(source + "\n")

    return msgpack_path, sources_path


def _discover_completed_sources(batch_dir: Path) -> Set[str]:
    """
    Resume quickly by reading only the small .sources files, not the full msgpack data.
    """
    completed: Set[str] = set()
    if not batch_dir.exists():
        return completed

    for source_file in sorted(batch_dir.glob("batch_*.sources")):
        with source_file.open("r", encoding="utf-8") as handle:
            for line in handle:
                source = line.strip()
                if source:
                    completed.add(source)

    return completed


def _next_batch_index(batch_dir: Path) -> int:
    """
    Continue numbering batch files correctly on resume.
    """
    if not batch_dir.exists():
        return 1

    max_index = 0
    for path in batch_dir.glob("batch_*.msgpack"):
        match = re.match(r"batch_(\d+)\.msgpack$", path.name)
        if match:
            max_index = max(max_index, int(match.group(1)))
    return max_index + 1


def _merge_batches(batch_dir: Path) -> Dict[str, Dict[str, object]]:
    """
    Merge all batch MessagePack files into one in-memory dictionary.

    This is the most memory-intensive step, but it happens only once at the end.
    """
    merged: Dict[str, Dict[str, object]] = {}

    batch_files = sorted(batch_dir.glob("batch_*.msgpack"))
    if not batch_files:
        return merged

    for batch_file in batch_files:
        batch_data = load_dict_msgpack(str(batch_file))
        merged.update(batch_data)
        del batch_data

    gc.collect()
    return merged


def _write_final_msgpack(path: Path, data: Dict[str, Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    save_dict_msgpack(data, str(path))


def _write_final_json(path: Path, data: Dict[str, Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    ordered = {key: data[key] for key in sorted(data)}
    with path.open("w", encoding="utf-8") as handle:
        json.dump(ordered, handle, indent=4, sort_keys=True)


def _record_topology_summary(record: Dict[str, object]) -> str:
    """
    One-line topology summary for progress output. For an "all" record this is
    one ``method=<net>`` item per node definition; otherwise the net.
    """
    topologies = record.get("topologies")
    if isinstance(topologies, dict):
        return " ".join(
            f"{method}={topologies.get(method, {}).get('topology')}"
            for method in ALL_METHODS
        )
    return str(record.get("topology"))


def _write_final_csv(path: Path, data: Dict[str, Dict[str, object]]) -> None:
    """
    Summary CSV only. Heavy CGD text fields are intentionally not written here.

    For "all" records each node definition gets its own columns
    (``<method>_topology`` and so on), one row per structure, so the table loads
    straight into a database.
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    is_all = any(
        isinstance(row.get("topologies"), dict) for row in data.values()
    )

    if is_all:
        fields = ("topology", "dimension", "td10", "topology_hash")
        fieldnames = ["record_id", "MOF name"]
        for method in ALL_METHODS:
            fieldnames.extend(f"{method}_{field}" for field in fields)
    else:
        fieldnames = [
            "record_id",
            "MOF name",
            "topology",
            "dimension",
            "td10",
            "topology_hash",
        ]

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()

        for key in sorted(data):
            row = data[key]
            out = {"record_id": key, "MOF name": row.get("mof_name")}
            if is_all:
                topologies = row.get("topologies") or {}
                for method in ALL_METHODS:
                    net = topologies.get(method, {})
                    for field in ("topology", "dimension", "td10", "topology_hash"):
                        out[f"{method}_{field}"] = net.get(field)
            else:
                for field in ("topology", "dimension", "td10", "topology_hash"):
                    out[field] = row.get(field)
            writer.writerow(out)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    batch_dir = Path(args.batch_dir)
    msgpack_path = Path(args.msgpack_out) if args.msgpack_out else None
    json_path = Path(args.json_out) if args.json_out else None
    csv_path = Path(args.csv_out) if args.csv_out else None

    if args.finalise_only:
        merged = _merge_batches(batch_dir)
        if not merged:
            raise SystemExit(f"No batch MessagePack files found in: {batch_dir}")

        if msgpack_path:
            _write_final_msgpack(msgpack_path, merged)
            print(f"Wrote MessagePack: {msgpack_path}")

        if json_path:
            _write_final_json(json_path, merged)
            print(f"Wrote JSON: {json_path}")

        if csv_path:
            _write_final_csv(csv_path, merged)
            print(f"Wrote CSV: {csv_path}")

        del merged
        gc.collect()
        return 0

    files = _collect_input_files(args.inputs, recursive=args.recursive)
    if not files:
        raise SystemExit("No input files were found.")

    completed_sources: Set[str] = set()
    if args.resume:
        completed_sources = _discover_completed_sources(batch_dir)

    next_batch_index = _next_batch_index(batch_dir) if args.resume else 1

    total_inputs = len(files)
    skipped = 0
    processed_new = 0

    batch_results: Dict[str, Dict[str, object]] = {}
    batch_sources: List[str] = []

    print(f"Found {total_inputs} candidate files.")
    if args.resume:
        print(f"Resume mode: {len(completed_sources)} previously flushed files will be skipped.")
    print(f"Starting batch index: {next_batch_index}")

    for filepath in files:
        resolved = str(filepath.resolve())

        if args.resume and resolved in completed_sources:
            skipped += 1
            continue

        key = _make_record_key(filepath)

        try:
            topo_data = _compute_topology(
                filepath,
                method=args.method,
                decimals=args.decimals,
                include_edge_centers=args.include_edge_centers,
                fallback_to_input_cgd=args.fallback_to_input_cgd,
            )
        except Exception:
            if args.method == "all":
                topo_data = {
                    "topologies": {m: dict(_ERROR_TOPOLOGY) for m in ALL_METHODS}
                }
            else:
                topo_data = dict(_ERROR_TOPOLOGY)

        record = {
            "mof_name": filepath.stem,
            "source": resolved,
            **topo_data,
        }

        batch_results[key] = record
        batch_sources.append(resolved)
        processed_new += 1

        print(f"{filepath}\t{_record_topology_summary(record)}")
        if args.verbose:
            print(json.dumps(record, indent=2, sort_keys=True))

        if len(batch_results) >= args.flush_every:
            msgpack_batch, sources_batch = _write_batch(
                batch_dir=batch_dir,
                batch_index=next_batch_index,
                batch_data=batch_results,
                batch_sources=batch_sources,
            )
            print(f"Flushed batch {next_batch_index}: {msgpack_batch}")
            print(f"Checkpoint file: {sources_batch}")

            completed_sources.update(batch_sources)
            batch_results.clear()
            batch_sources.clear()
            next_batch_index += 1
            gc.collect()

        del topo_data
        del record

    if batch_results:
        msgpack_batch, sources_batch = _write_batch(
            batch_dir=batch_dir,
            batch_index=next_batch_index,
            batch_data=batch_results,
            batch_sources=batch_sources,
        )
        print(f"Flushed final batch {next_batch_index}: {msgpack_batch}")
        print(f"Checkpoint file: {sources_batch}")

        completed_sources.update(batch_sources)
        batch_results.clear()
        batch_sources.clear()
        gc.collect()

    print(f"Processing complete.")
    print(f"Newly processed this run: {processed_new}")
    print(f"Skipped from checkpoint: {skipped}")

    merged = _merge_batches(batch_dir)
    if not merged:
        raise SystemExit(f"No batch MessagePack files found in: {batch_dir}")

    if msgpack_path:
        _write_final_msgpack(msgpack_path, merged)
        print(f"Wrote MessagePack: {msgpack_path}")

    if json_path:
        _write_final_json(json_path, merged)
        print(f"Wrote JSON: {json_path}")

    if csv_path:
        _write_final_csv(csv_path, merged)
        print(f"Wrote CSV: {csv_path}")

    print(f"Done. Total merged records: {len(merged)}")
    print(f"Intermediate batches remain in: {batch_dir}")

    del merged
    gc.collect()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
