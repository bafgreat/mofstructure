#!/usr/bin/env python3
"""
mofstructure_systre CLI

Works out of the box for:
  - .cgd files (run Systre directly)
  - structure files readable by ASE (e.g., .cif) -> generate CGD then run Systre
  - folders containing mix of files (runs per file)

Folder output:
  - Optionally write CSV and/or JSON:
      --csv results.csv
      --json results.json
    where the key is the file basename (stem) and the value is topology.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, Tuple

from mofstructure.systre import identify_topology, identify_topology_batch


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Run Systre (Gavrog) on CGD files or structures (auto-generate CGD when needed)."
    )
    p.add_argument("input", help="Path to a .cgd file, a structure file (e.g. .cif), or a folder.")

    # systre/java
    p.add_argument("--java", default="java", help="Java executable (default: java).")
    p.add_argument("--timeout", type=int, default=30, help="Systre timeout in seconds (default: 30).")
    p.add_argument("--xmx", default="1024m", help="Java max heap, e.g. 1024m (default: 1024m).")
    p.add_argument("--keep-tmp", action="store_true", help="Keep temporary CGD files (debug).")

    # folder options
    p.add_argument("--recursive", action="store_true", default=True, help="Recurse into folders (default: True).")
    p.add_argument("--no-recursive", action="store_false", dest="recursive", help="Do not recurse into folders.")

    # CGD generation options (used when input is NOT .cgd)
    p.add_argument(
        "--method",
        choices=["sbus", "all_node", "ligand_cluster"],
        default="all_node",
        help="CGD method for structure inputs (default: all_node).",
    )
    p.add_argument("--name", default="net", help="CGD graph ID name (default: net).")
    p.add_argument("--top-k-regions", type=int, default=1, help="For all_node: pick top K regions (default: 1).")
    p.add_argument("--connect-mode", choices=["clique", "chain"], default="clique", help="For all_node contraction.")
    p.add_argument(
        "--dedup-mode",
        choices=["none", "shift", "topological"],
        default="shift",
        help="Dedup mode for CGD generation (default: shift).",
    )
    p.add_argument("--no-dedup", action="store_true", help="Disable edge dedup during CGD generation.")

    # batch filtering
    p.add_argument(
        "--patterns",
        nargs="*",
        default=None,
        help="Optional list of file extensions to consider in folders (e.g. .cgd .cif .vasp). "
             "Default: auto (common structure formats + .cgd).",
    )

    # outputs (folder mode; also allowed in single-file mode)
    p.add_argument("--csv", dest="csv_out", default='Topology.csv', help="Write results to CSV file (basename, topology).")
    p.add_argument("--json", dest="json_out", default='Topology.json', help="Write results to JSON file ({basename: topology}).")

    # output verbosity
    p.add_argument("--verbose", action="store_true", help="Print Systre stdout/stderr.")
    return p


def _write_csv(path: Path, mapping: Dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["MOF name", "topology"])
        for k in sorted(mapping.keys()):
            w.writerow([k, mapping[k]])


def _write_json(path: Path, mapping: Dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    # stable ordering for nicer diffs
    ordered = {k: mapping[k] for k in sorted(mapping.keys())}
    path.write_text(json.dumps(ordered, indent=4, sort_keys=True), encoding="utf-8")


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    inp = Path(args.input)

    common = dict(
        java=args.java,
        timeout_s=args.timeout,
        xmx=args.xmx,
        keep_tmp=args.keep_tmp,
    )

    gen = dict(
        method=args.method,
        name=args.name,
        top_k_regions=args.top_k_regions,
        connect_mode=args.connect_mode,
        dedup_parallel_edges=(not args.no_dedup),
        dedup_mode=args.dedup_mode,
    )

    patterns: Tuple[str, ...] = tuple(args.patterns) if args.patterns else (
        ".cgd", ".cif", ".vasp", ".poscar", ".xyz", ".pdb", ".mol", ".sdf"
    )

    if not inp.exists():
        raise SystemExit(f"Input not found: {inp}")

    # Folder batch mode
    if inp.is_dir():
        resmap = identify_topology_batch(
            inp,
            patterns=patterns,
            recursive=args.recursive,
            **gen,
            **common,
        )
      
        topo_map: Dict[str, str] = {}

        # handle possible basename collisions (same stem in different subfolders)
        # We keep "stem" as key, but if collision happens, we suffix with __2, __3...
        counts: Dict[str, int] = {}

        for fp in sorted(resmap.keys()):
            res = resmap[fp]
            stem = Path(fp).stem
            if stem in topo_map:
                counts[stem] = counts.get(stem, 1) + 1
                stem_key = f"{stem}__{counts[stem]}"
            else:
                stem_key = stem
                counts[stem] = 1

            topo_map[stem_key] = res.topology

            print(f"{fp}\t{res.topology}")
            if args.verbose and (res.stdout or res.stderr):
                print("---- STDOUT ----")
                print(res.stdout)
                print("---- STDERR ----")
                print(res.stderr)
                print("----------------")

        # write outputs if requested
        if args.csv_out:
            _write_csv(Path(args.csv_out), topo_map)
        if args.json_out:
            _write_json(Path(args.json_out), topo_map)

        return 0

    # Single file mode
    if inp.suffix.lower() == ".cgd":
        res = identify_topology(inp, input_is_cgd=True, **common)
    else:
        res = identify_topology(inp, **gen, **common)

    print(res.topology)

    # if user requested csv/json for a single file, write a single-entry mapping
    if args.csv_out or args.json_out:
        topo_map = {inp.stem: res.topology}
        if args.csv_out:
            _write_csv(Path(args.csv_out), topo_map)
        if args.json_out:
            _write_json(Path(args.json_out), topo_map)

    if args.verbose and (res.stdout or res.stderr):
        print("---- STDOUT ----")
        print(res.stdout)
        print("---- STDERR ----")
        print(res.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
