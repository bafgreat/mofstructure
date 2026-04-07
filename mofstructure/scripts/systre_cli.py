#!/usr/bin/env python3
"""
mofstructure_topology CLI

This CLI uses `MOFstructure.get_topology()` from `mofstructure.structure`.

Supported input:
  - one file
  - multiple files
  - one or more folders

Behavior:
  - if an input is a file, it is processed directly
  - if an input is a folder, all ASE-readable files inside it are found and processed
  - results can be written to CSV and/or JSON

Notes
-----
This script delegates topology detection to:

    MOFstructure(filename=...).get_topology(...)

Therefore, it automatically uses the topology workflow already implemented
in `mofstructure.structure` and `mofstructure.systre`.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from ase.io import read
from mofstructure.structure import MOFstructure


def build_parser() -> argparse.ArgumentParser:
    """
    Build the command-line parser.

    **returns:**
        argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser(
        description="Compute topology using MOFstructure.get_topology() for files or folders."
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="One or more files and/or folders.",
    )
    parser.add_argument(
        "--method",
        choices=["sbus", "all_node", "ligand_cluster"],
        default="all_node",
        help="Topology method passed to MOFstructure.get_topology().",
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
        help="Fallback to input CGD text if no relaxed component is parsed.",
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
        "--csv",
        dest="csv_out",
        default="Topology.csv",
        help="Write results to CSV file.",
    )
    parser.add_argument(
        "--json",
        dest="json_out",
        default="Topology.json",
        help="Write results to JSON file.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print full topology dictionaries.",
    )
    return parser


def _is_ase_readable(path: Path) -> bool:
    """
    Check whether ASE can read a file.

    **parameters:**
        path: pathlib.Path
            Input file path.

    **returns:**
        bool
    """
    try:
        read(path)
        return True
    except Exception:
        return False


def _collect_input_files(inputs: Sequence[str], recursive: bool = True) -> List[Path]:
    """
    Collect all readable files from input paths.

    Rules:
      - files are added directly
      - directories are searched for ASE-readable files

    **parameters:**
        inputs: sequence
            Input file and/or folder paths.

        recursive: bool
            If True, recurse into folders.

    **returns:**
        list
            List of unique readable file paths.
    """
    collected: List[Path] = []
    seen = set()

    for item in inputs:
        path = Path(item)

        if not path.exists():
            print(f"Warning: input not found: {path}")
            continue

        if path.is_file():
            if _is_ase_readable(path):
                resolved = str(path.resolve())
                if resolved not in seen:
                    collected.append(path)
                    seen.add(resolved)
            else:
                print(f"Warning: ASE could not read file: {path}")
            continue

        if path.is_dir():
            iterator = path.rglob("*") if recursive else path.glob("*")
            for subpath in iterator:
                if not subpath.is_file():
                    continue
                if not _is_ase_readable(subpath):
                    continue
                resolved = str(subpath.resolve())
                if resolved not in seen:
                    collected.append(subpath)
                    seen.add(resolved)

    return collected


def _compute_topology(
    path: Path,
    *,
    method: str,
    decimals: int,
    include_edge_centers: bool,
    fallback_to_input_cgd: bool,
) -> Dict[str, object]:
    """
    Compute topology data for one structure file.

    **parameters:**
        path: pathlib.Path
            Input structure path.

        method: str
            Topology method.

        decimals: int
            Number of decimal places used in topology hashing.

        include_edge_centers: bool
            If True, include edge-center comments in the CRYSTAL2 text.

        fallback_to_input_cgd: bool
            If True, fallback to input CGD text when needed.

    **returns:**
        python dictionary
            Topology data returned by `MOFstructure.get_topology()`.
    """
    structure = MOFstructure(filename=str(path))
    data = structure.get_topology(
        method=method,
        decimals=decimals,
        include_edge_centers=include_edge_centers,
        fallback_to_input_cgd=fallback_to_input_cgd,
    )
    return data


def _write_csv(path: Path, mapping: Dict[str, Dict[str, object]]) -> None:
    """
    Write topology results to CSV.

    **parameters:**
        path: pathlib.Path
            Output CSV path.

        mapping: python dictionary
            Mapping:
                basename -> topology data dictionary
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "MOF name",
        "topology",
        "dimension",
        "td10",
        "topology_hash",
        "cgd_crystal2text",
    ]

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()

        for key in sorted(mapping):
            row = mapping[key]
            writer.writerow(
                {
                    "MOF name": key,
                    "topology": row.get("topology"),
                    "dimension": row.get("dimension"),
                    "td10": row.get("td10"),
                    "topology_hash": row.get("topology_hash"),
                    "cgd_crystal2text": row.get("cgd_crystal2text"),
                }
            )


def _write_json(path: Path, mapping: Dict[str, Dict[str, object]]) -> None:
    """
    Write topology results to JSON.

    **parameters:**
        path: pathlib.Path
            Output JSON path.

        mapping: python dictionary
            Mapping:
                basename -> topology data dictionary
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    ordered = {key: mapping[key] for key in sorted(mapping)}
    path.write_text(json.dumps(ordered, indent=4, sort_keys=True), encoding="utf-8")


def main(argv: Optional[Sequence[str]] = None) -> int:
    """
    Command-line entry point.

    **parameters:**
        argv: list, optional
            Command-line arguments.

    **returns:**
        int
            Exit status code.
    """
    args = build_parser().parse_args(argv)

    files = _collect_input_files(args.inputs, recursive=args.recursive)
    if not files:
        raise SystemExit("No ASE-readable files were found.")

    results: Dict[str, Dict[str, object]] = {}
    counts: Dict[str, int] = {}

    for filepath in files:
        try:
            topo_data = _compute_topology(
                filepath,
                method=args.method,
                decimals=args.decimals,
                include_edge_centers=args.include_edge_centers,
                fallback_to_input_cgd=args.fallback_to_input_cgd,
            )
        except Exception as exc:
            topo_data = {
                "topology": "ERROR",
                "dimension": None,
                "td10": None,
                "topology_hash": None,
                "cgd_crystal2text": None,
                "error": str(exc),
            }

        stem = filepath.stem
        if stem in results:
            counts[stem] = counts.get(stem, 1) + 1
            key = f"{stem}__{counts[stem]}"
        else:
            counts[stem] = 1
            key = stem

        results[key] = topo_data

        print(f"{filepath}\t{topo_data.get('topology')}")
        if args.verbose:
            print(json.dumps(topo_data, indent=2, sort_keys=True))

    if args.csv_out:
        _write_csv(Path(args.csv_out), results)

    if args.json_out:
        _write_json(Path(args.json_out), results)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())