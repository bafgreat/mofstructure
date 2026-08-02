#!/usr/bin/env python3
from __future__ import annotations
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"

import argparse
from pathlib import Path
from typing import Optional, Sequence

from mofstructure.systre import SystreTopology


def build_argparser() -> argparse.ArgumentParser:
    """
    Build the command-line argument parser.

    **returns:**
        argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser(
        description=(
            "Generate a PERIODIC_GRAPH CGD with mofstructure, run Systre, "
            "save stdout/stderr, and write a CRYSTAL2 CGD from the best relaxed component."
        )
    )
    parser.add_argument("input_structure", help="Input structure file, e.g. CIF.")
    parser.add_argument("--name", default="net", help="Graph name used for CGD generation.")
    parser.add_argument(
        "--method",
        default="sbus",
        help="Topology extraction method passed to SystreTopology.",
    )
    parser.add_argument(
        "--systre-log",
        default="systre.log",
        help="Path to write Systre stdout.",
    )
    parser.add_argument(
        "--systre-err",
        default="systre.stderr",
        help="Path to write Systre stderr.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="systre_crystal2.cgd",
        help="Output CRYSTAL2 CGD file.",
    )
    parser.add_argument(
        "--keep-periodic-graph",
        default=None,
        help="Optional path to save the PERIODIC_GRAPH CGD used as input to Systre.",
    )
    parser.add_argument(
        "--include-edge-centers",
        action="store_true",
        help="Include commented edge-center lines in the CRYSTAL2 output.",
    )
    parser.add_argument(
        "--fallback-to-input-cgd",
        action="store_true",
        help="Fallback to the input PERIODIC_GRAPH CGD if relaxed output is unavailable.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print where the Systre log and error files were written.",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """
    Run the Systre workflow from the command line.

    **parameters:**
        argv: list, optional
            Command-line arguments.

    **returns:**
        int
            Exit status code.
    """
    args = build_argparser().parse_args(argv)

    input_path = Path(args.input_structure)
    if not input_path.exists():
        raise FileNotFoundError(input_path)

    runner = SystreTopology(
        input_path,
        method=args.method,
        name=args.name,
        keep_tmp=False,
    )

    result = runner.identify()

    Path(args.systre_log).write_text(result.stdout or "", encoding="utf-8")
    Path(args.systre_err).write_text(result.stderr or "", encoding="utf-8")

    if args.verbose:
        print(f"Systre log written to {args.systre_log}")
        print(f"Systre errors written to {args.systre_err}")

    if args.keep_periodic_graph is not None:
        cgd_text = runner._input_cgd_text() if hasattr(runner, "_input_cgd_text") else None
        if cgd_text:
            Path(args.keep_periodic_graph).write_text(cgd_text, encoding="utf-8")

    best = runner.best_component()

    if best is not None:
        runner.write_crystal2(
            args.output,
            component_index=best.component_index,
            include_edge_centers=args.include_edge_centers,
            fallback_to_input_cgd=args.fallback_to_input_cgd,
        )
        print(f"Topology: {result.topology}")
        print(f"RCSR name: {best.rcsr_name}")
        print(f"Dimension: {best.dimension}")
        print(f"TD10: {best.td10}")
        print(f"Topology hash: {best.topology_hash()}")
        print(f"Wrote CRYSTAL2: {args.output}")
    else:
        crystal2_text = runner.crystal2_text(
            include_edge_centers=args.include_edge_centers,
            fallback_to_input_cgd=args.fallback_to_input_cgd,
        )
        if crystal2_text is None:
            raise RuntimeError(
                "Systre completed, but no relaxed component could be parsed and no fallback CRYSTAL2 text was available."
            )

        Path(args.output).write_text(crystal2_text, encoding="utf-8")
        print(f"Topology: {result.topology}")
        print("No relaxed component parsed; wrote fallback CRYSTAL2 text.")
        print(f"Wrote CRYSTAL2: {args.output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())