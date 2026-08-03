#!/usr/bin/env python3
"""Command-line interface for ligand--cluster MOF fingerprints."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set

from mofstructure.filetyper import convert_numpy_types
from mofstructure.structure import MOFstructure


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compute ligand--cluster fingerprints for MOF files or folders."
    )
    parser.add_argument("inputs", nargs="+", help="One or more structure files or folders.")
    parser.add_argument(
        "--recursive", action="store_true", default=True,
        help="Recurse into folders (default: true).",
    )
    parser.add_argument(
        "--no-recursive", action="store_false", dest="recursive",
        help="Do not recurse into folders.",
    )
    parser.add_argument(
        "--json", dest="json_out", default="Fingerprints.json",
        help="Output JSON path (default: Fingerprints.json).",
    )
    parser.add_argument(
        "--csv", dest="csv_out", default="Fingerprints.csv",
        help="Output CSV path (default: Fingerprints.csv).",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    return parser


def _collect_input_files(inputs: Sequence[str], recursive: bool = True) -> List[Path]:
    files: List[Path] = []
    seen: Set[str] = set()
    for item in inputs:
        path = Path(item)
        if not path.exists():
            print(f"Warning: input not found: {path}")
            continue
        candidates = [path] if path.is_file() else (
            path.rglob("*") if recursive else path.glob("*")
        )
        for candidate in candidates:
            if not candidate.is_file():
                continue
            resolved = str(candidate.resolve())
            if resolved not in seen:
                files.append(candidate)
                seen.add(resolved)
    return files


def _record_key(path: Path) -> str:
    stem = re.sub(r"[^A-Za-z0-9_.-]", "_", re.sub(r"\s+", "_", path.stem.strip()))
    digest = hashlib.sha1(str(path.resolve()).encode("utf-8")).hexdigest()[:12]
    return f"{stem or 'structure'}__{digest}"


def _compute_fingerprint(path: Path) -> Dict[str, object]:
    return convert_numpy_types(
        MOFstructure(filename=str(path)).get_ligand_cluster_fingerprint()
    )


def _write_json(path: Path, records: Dict[str, Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump({key: records[key] for key in sorted(records)}, handle,
                  indent=4, sort_keys=True)


def _write_csv(path: Path, records: Dict[str, Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["record_id", "MOF name", "source", "fingerprint_hash", "cluster_units"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for key in sorted(records):
            record = records[key]
            writer.writerow({
                "record_id": key,
                "MOF name": record.get("mof_name"),
                "source": record.get("source"),
                "fingerprint_hash": record.get("fingerprint_hash"),
                "cluster_units": record.get("cluster_units"),
            })


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    files = _collect_input_files(args.inputs, recursive=args.recursive)
    if not files:
        raise SystemExit("No input files were found.")

    records: Dict[str, Dict[str, object]] = {}
    failures = 0
    for path in files:
        try:
            fingerprint = _compute_fingerprint(path)
        except Exception as error:
            failures += 1
            print(f"Error: {path}: {error}")
            continue
        record = {
            "mof_name": path.stem,
            "source": str(path.resolve()),
            **fingerprint,
        }
        records[_record_key(path)] = record
        print(f"{path}\t{record['fingerprint_hash']}")
        if args.verbose:
            print(json.dumps(record, indent=2, sort_keys=True))

    if not records:
        raise SystemExit("Fingerprinting failed for every input file.")
    _write_json(Path(args.json_out), records)
    _write_csv(Path(args.csv_out), records)
    print(f"Wrote JSON: {args.json_out}")
    print(f"Wrote CSV: {args.csv_out}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
