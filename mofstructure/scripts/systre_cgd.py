#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import tempfile
from contextlib import contextmanager
from dataclasses import dataclass
from importlib.resources import as_file, files as ir_files
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple
from mofstructure import generate_cgd


# =========================
# Java / Systre
# =========================

def find_java(java: str = "java") -> str:
    """Find a working java executable."""
    # explicit path
    if os.path.sep in java or (os.path.altsep and os.path.altsep in java):
        p = Path(os.path.expandvars(os.path.expanduser(java)))
        if p.exists():
            return str(p)
        raise RuntimeError(f"Java executable not found: {p}")

    exe = shutil.which(java)
    if exe:
        return exe

    java_home = os.environ.get("JAVA_HOME")
    if java_home:
        candidate = Path(java_home) / "bin" / ("java.exe" if os.name == "nt" else "java")
        if candidate.exists():
            return str(candidate)

    raise RuntimeError(
        "Java not found. Install Java and ensure `java` is on PATH, "
        "or pass --java /path/to/java, or set JAVA_HOME."
    )


def run_systre_on_cgd(
    cgd_path: Path,
    *,
    jar: Path,
    arc: Path,
    java: str = "java",
    xmx: str = "1024m",
    timeout_s: int = 60,
) -> subprocess.CompletedProcess:
    java_exe = find_java(java)
    cmd = [
        java_exe,
        f"-Xmx{xmx}",
        "-cp",
        str(jar),
        "org.gavrog.apps.systre.SystreCmdline",
        str(arc),
        str(cgd_path),
    ]
    return subprocess.run(
        cmd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=timeout_s,
        check=False,
    )


# =========================
# Installed-package resources
# =========================

@contextmanager
def systre_resources(jar_arg: str, arc_arg: str):
    """
    Yield (jar_path, arc_path) as real filesystem Paths.

    Behavior:
    - If user passes an existing path, use it.
    - Else, load from installed `mofstructure` package resources:
        JAR: mofstructure/bin/Systre-19.6.0.jar
        ARC: mofstructure/db/RCSRnets-2019-06-01.arc
    """

    def expand(s: str) -> Path:
        return Path(os.path.expandvars(os.path.expanduser(s)))

    def strip_pkg_prefix(s: str) -> str:
        s = os.path.expandvars(os.path.expanduser(s))
        prefix = "mofstructure" + os.sep
        return s[len(prefix):] if s.startswith(prefix) else s

    @contextmanager
    def resource_or_path(arg: str, fallback_rel: str):
        p = expand(arg)
        if p.exists():
            yield p.resolve()
            return

        rel = strip_pkg_prefix(arg)
        cand = ir_files("mofstructure").joinpath(rel)
        if not cand.exists():
            cand = ir_files("mofstructure").joinpath(fallback_rel)

        if not cand.exists():
            raise FileNotFoundError(
                "Systre resource not found.\n"
                f"  Tried filesystem path: {p}\n"
                f"  Tried package resource: mofstructure/{rel}\n"
                f"  Tried fallback resource: mofstructure/{fallback_rel}\n\n"
                "Fix:\n"
                "  - Ensure these files are included as package-data in your wheel/sdist, OR\n"
                "  - Provide absolute paths via --jar and --arc.\n"
            )

        with as_file(cand) as real:
            yield Path(real).resolve()

    with resource_or_path(jar_arg, "bin/Systre-19.6.0.jar") as jar:
        with resource_or_path(arc_arg, "db/RCSRnets-2019-06-01.arc") as arc:
            yield jar, arc


# =========================
# Parse relaxed output
# =========================

@dataclass
class RelaxedComponent:
    component_index: int
    dimension: int
    cell: Tuple[float, float, float, float, float, float]  # a b c alpha beta gamma
    nodes: Dict[str, Tuple[float, float, float]]
    edges: List[Tuple[Tuple[float, float, float], Tuple[float, float, float]]]
    edge_centers: List[Tuple[float, float, float]]
    coordination_number: Dict[str, int]
    given_space_group: Optional[str] = None
    ideal_space_group: Optional[str] = None
    rcsr_name: Optional[str] = None
    td10: Optional[int] = None

def _split_component_blocks(stdout: str) -> List[Tuple[int, str]]:
    if "Processing component " not in stdout:
        return [(1, stdout)]
    parts = re.split(r"\n\s*Processing component\s+(\d+):\s*\n", stdout)
    if len(parts) < 3:
        return [(1, stdout)]
    pre = parts[0]
    out: List[Tuple[int, str]] = []
    for i in range(1, len(parts) - 1, 2):
        idx = int(parts[i])
        block = parts[i + 1]
        out.append((idx, pre + "\nProcessing component %d:\n" % idx + block))
    return out


def _split_structure_blocks(stdout: str) -> List[Tuple[int, str]]:
    """
    Split Systre stdout into per-structure blocks.

    Supports:
      - "Structure #1 - "net". ... Finished structure #1 - "net"."
      - fallback: whole stdout as one block
    """
    m = list(re.finditer(r"(?m)^\s*Structure\s+#(\d+)\s+-\s+\".*?\"\.\s*$", stdout))
    if not m:
        return [(1, stdout)]

    blocks: List[Tuple[int, str]] = []
    for i, mi in enumerate(m):
        idx = int(mi.group(1))
        start = mi.start()
        end = m[i + 1].start() if i + 1 < len(m) else len(stdout)
        blocks.append((idx, stdout[start:end]))
    return blocks

def _parse_td10(block: str) -> Optional[int]:
    m = re.search(r"(?m)^\s*TD10\s*=\s*(\d+)\s*$", block)
    return int(m.group(1)) if m else None


def _parse_rcsr_name(block: str) -> Optional[str]:
    """
    Prefer the name under 'Structure was identified with RCSR symbol:'
    Example:
        Structure was identified with RCSR symbol:
            Name:        pcu
    """
    m = re.search(
        r"Structure\s+was\s+identified\s+with\s+RCSR\s+symbol:\s*(?:\r?\n)+\s*Name:\s*([A-Za-z0-9_\-]+)",
        block,
        flags=re.IGNORECASE,
    )
    if m:
        return m.group(1)

    # Fallback: sometimes only the "found in archive" block is present
    m2 = re.search(
        r"Structure\s+was\s+found\s+in\s+archive.*?(?:\r?\n)+\s*Name:\s*([A-Za-z0-9_\-]+)",
        block,
        flags=re.IGNORECASE | re.DOTALL,
    )
    return m2.group(1) if m2 else None

def _parse_dimension(block: str) -> Optional[int]:
    # "Input structure described as 3-periodic."
    m = re.search(r"described\s+as\s+(\d+)\s*-\s*periodic", block, flags=re.IGNORECASE)
    if m:
        return int(m.group(1))

    # "Input structure described as 2 periodic."
    m = re.search(r"described\s+as\s+(\d+)\s+periodic", block, flags=re.IGNORECASE)
    if m:
        return int(m.group(1))

    # old style "dimension = 3"
    m = re.search(r"dimension\s*=\s*(\d+)", block, flags=re.IGNORECASE)
    if m:
        return int(m.group(1))

    return None


def _parse_relaxed_cell(block: str) -> Optional[Tuple[float, float, float, float, float, float]]:
    """
    Return (a,b,c,alpha,beta,gamma).
    For 2D, set c=1.0, alpha=beta=90, gamma=parsed.
    """

    # Your multi-line 3D format:
    # Relaxed cell parameters:
    #   a = ..., b = ..., c = ...
    #   alpha = ..., beta = ..., gamma = ...
    m3 = re.search(
        r"Relaxed cell parameters:\s*(?:\r?\n)\s*"
        r"a\s*=\s*([0-9.\-Ee+]+)\s*,\s*b\s*=\s*([0-9.\-Ee+]+)\s*,\s*c\s*=\s*([0-9.\-Ee+]+)\s*(?:\r?\n)\s*"
        r"alpha\s*=\s*([0-9.\-Ee+]+)\s*,\s*beta\s*=\s*([0-9.\-Ee+]+)\s*,\s*gamma\s*=\s*([0-9.\-Ee+]+)",
        block,
        flags=re.IGNORECASE,
    )
    if m3:
        return (
            float(m3.group(1)),
            float(m3.group(2)),
            float(m3.group(3)),
            float(m3.group(4)),
            float(m3.group(5)),
            float(m3.group(6)),
        )

    # Common 2D format (some Systre builds):
    # Relaxed cell parameters:
    #   a = ..., b = ..., gamma = ...
    m2 = re.search(
        r"Relaxed cell parameters:\s*(?:\r?\n)\s*"
        r"a\s*=\s*([0-9.\-Ee+]+)\s*,\s*b\s*=\s*([0-9.\-Ee+]+)\s*,\s*gamma\s*=\s*([0-9.\-Ee+]+)",
        block,
        flags=re.IGNORECASE,
    )
    if m2:
        a = float(m2.group(1))
        b = float(m2.group(2))
        gamma = float(m2.group(3))
        return (a, b, 1.0, 90.0, 90.0, gamma)

    # Older single-line formats (keep as fallback)
    m1 = re.search(
        r"Relaxed cell parameters:\s*(?:\r?\n|\s+)"
        r"a\s*=\s*([0-9.\-Ee+]+)\s*,\s*b\s*=\s*([0-9.\-Ee+]+)\s*,\s*c\s*=\s*([0-9.\-Ee+]+)\s*"
        r"(?:,|\s+)\s*alpha\s*=\s*([0-9.\-Ee+]+)\s*,\s*beta\s*=\s*([0-9.\-Ee+]+)\s*,\s*gamma\s*=\s*([0-9.\-Ee+]+)",
        block,
        flags=re.IGNORECASE,
    )
    if m1:
        return (
            float(m1.group(1)),
            float(m1.group(2)),
            float(m1.group(3)),
            float(m1.group(4)),
            float(m1.group(5)),
            float(m1.group(6)),
        )

    return None

def _parse_space_groups(block: str) -> Tuple[Optional[str], Optional[str]]:
    # "Given space group is P1."
    m_given = re.search(r"Given\s+space\s+group\s+is\s+([A-Za-z0-9\-_/]+)\s*\.", block, flags=re.IGNORECASE)
    given = m_given.group(1) if m_given else None

    # "Ideal space group is Pm-3m."
    m_ideal = re.search(r"Ideal\s+space\s+group\s+is\s+([A-Za-z0-9\-_/]+)\s*\.", block, flags=re.IGNORECASE)
    ideal = m_ideal.group(1) if m_ideal else None

    return given, ideal


def _parse_edge_centers(block: str) -> List[Tuple[float, float, float]]:
    """
    Parse:
      Edge centers:
         x y z
         x y z
    Stops at next header.
    """
    m = re.search(
        r"Edge\s+centers:\s*(.*?)(?:\r?\n)\s*(?:Edge statistics:|Angle statistics:|Shortest non-bonded distance|Degrees of freedom:|Finished structure|\Z)",
        block,
        flags=re.DOTALL | re.IGNORECASE,
    )
    if not m:
        return []

    centers: List[Tuple[float, float, float]] = []
    for line in m.group(1).strip().splitlines():
        line = line.strip()
        if not line:
            continue

        # 3D center: x y z
        m3 = re.match(r"([0-9.\-Ee+]+)\s+([0-9.\-Ee+]+)\s+([0-9.\-Ee+]+)", line)
        if m3:
            centers.append((float(m3.group(1)), float(m3.group(2)), float(m3.group(3))))
            continue

        # 2D center: x y  (store z=0)
        m2 = re.match(r"([0-9.\-Ee+]+)\s+([0-9.\-Ee+]+)", line)
        if m2:
            centers.append((float(m2.group(1)), float(m2.group(2)), 0.0))

    return centers

def _parse_relaxed_positions(block: str) -> Optional[Dict[str, Tuple[float, float, float]]]:
    m_pos = re.search(r"Relaxed positions:\s*(.*?)(?:\r?\n)\s*Edges:", block, flags=re.DOTALL | re.IGNORECASE)
    if not m_pos:
        return None

    nodes: Dict[str, Tuple[float, float, float]] = {}
    for line in m_pos.group(1).strip().splitlines():
        line = line.strip()

        # 3D: Node 1: x y z
        m3 = re.match(r"Node\s+(\S+):\s+([0-9.\-Ee+]+)\s+([0-9.\-Ee+]+)\s+([0-9.\-Ee+]+)", line)
        if m3:
            nodes[m3.group(1)] = (float(m3.group(2)), float(m3.group(3)), float(m3.group(4)))
            continue

        # 2D: Node 1: x y
        m2 = re.match(r"Node\s+(\S+):\s+([0-9.\-Ee+]+)\s+([0-9.\-Ee+]+)", line)
        if m2:
            nodes[m2.group(1)] = (float(m2.group(2)), float(m2.group(3)), 0.0)

    return nodes

def _parse_coordination_numbers(block: str) -> Dict[str, int]:
    """
    Parse:
      Coordination sequences:
         Node 1:    6 18 38 ...
    Return: {"1": 6, ...}
    """
    m = re.search(
        r"Coordination\s+sequences:\s*(.*?)(?:\r?\n\s*\r?\n|\r?\n\s*TD10\s*=|\Z)",
        block,
        flags=re.DOTALL | re.IGNORECASE,
    )
    if not m:
        return {}

    out: Dict[str, int] = {}
    for line in m.group(1).splitlines():
        line = line.strip()
        # Node 1:    6 18 38 ...
        mm = re.match(r"Node\s+(\S+):\s+(\d+)\b", line)
        if mm:
            out[mm.group(1)] = int(mm.group(2))
    return out

def _parse_edges(block: str) -> Optional[List[Tuple[Tuple[float, float, float], Tuple[float, float, float]]]]:
    # Try to stop at "Edge centers:" if present; else other common section headers
    m_edges = re.search(
        r"Edges:\s*(.*?)(?:\r?\n)\s*(?:Edge centers:|Edge statistics:|Angle statistics:|Shortest non-bonded distance|Degrees of freedom:|Finished structure|\Z)",
        block,
        flags=re.DOTALL | re.IGNORECASE,
    )
    if not m_edges:
        return None

    edges: List[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = []
    for line in m_edges.group(1).strip().splitlines():
        line = line.strip()

        # 3D edge
        m3 = re.match(
            r"([0-9.\-Ee+]+)\s+([0-9.\-Ee+]+)\s+([0-9.\-Ee+]+)\s+<->\s+([0-9.\-Ee+]+)\s+([0-9.\-Ee+]+)\s+([0-9.\-Ee+]+)",
            line,
        )
        if m3:
            p = (float(m3.group(1)), float(m3.group(2)), float(m3.group(3)))
            q = (float(m3.group(4)), float(m3.group(5)), float(m3.group(6)))
            edges.append((p, q))
            continue

        # 2D edge
        m2 = re.match(
            r"([0-9.\-Ee+]+)\s+([0-9.\-Ee+]+)\s+<->\s+([0-9.\-Ee+]+)\s+([0-9.\-Ee+]+)",
            line,
        )
        if m2:
            p = (float(m2.group(1)), float(m2.group(2)), 0.0)
            q = (float(m2.group(3)), float(m2.group(4)), 0.0)
            edges.append((p, q))

    return edges


def parse_relaxed_component(block: str, component_index: int) -> Optional[RelaxedComponent]:
    dim = _parse_dimension(block) or 0

    cell = _parse_relaxed_cell(block)
    if not cell:
        return None

    nodes = _parse_relaxed_positions(block)
    if not nodes:
        return None

    edges = _parse_edges(block)
    if edges is None:
        return None

    edge_centers = _parse_edge_centers(block)
    given_sg, ideal_sg = _parse_space_groups(block)

    rcsr_name = _parse_rcsr_name(block)
    td10 = _parse_td10(block)
    coord_nums = _parse_coordination_numbers(block)

    return RelaxedComponent(
        component_index=component_index,
        dimension=dim,
        cell=cell,
        nodes=nodes,
        edges=edges,
        edge_centers=edge_centers,
        coordination_number=coord_nums,
        given_space_group=given_sg,
        ideal_space_group=ideal_sg,
        rcsr_name=rcsr_name,
        td10=td10
        )

def _sort_key(k: str):
    try:
        return (0, int(k))
    except ValueError:
        return (1, k)

# def write_crystal2_from_relaxed(out_path: Path, *, name: str, comp: RelaxedComponent) -> None:
#     """Write a CRYSTAL2 CGD from a parsed relaxed component."""
#     coord_to_id: Dict[Tuple[float, float, float], int] = {}
#     node_coords: List[Tuple[float, float, float]] = []
#     edge_centers: List[Tuple[float, float, float]] = comp.edge_centers
#     node_items = sorted(comp.nodes.items(), key=lambda kv: _sort_key(kv[0]))


#     for xyz in comp.nodes.values():
#         if xyz not in coord_to_id:
#             coord_to_id[xyz] = len(node_coords) + 1
#             node_coords.append(xyz)

#     for p, q in comp.edges:
#         for xyz in (p, q):
#             if xyz not in coord_to_id:
#                 coord_to_id[xyz] = len(node_coords) + 1
#                 node_coords.append(xyz)


#     seen = set()
#     clean_edges: List[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = []
#     for p, q in comp.edges:
#         key = (p, q) if p <= q else (q, p)
#         if key in seen:
#             continue
#         seen.add(key)
#         clean_edges.append((p, q))

#     deg = [0] * len(node_coords)
#     for p, q in clean_edges:
#         deg[coord_to_id[p] - 1] += 1
#         deg[coord_to_id[q] - 1] += 1

#     a, b, c, alpha, beta, gamma = comp.cell

#     lines: List[str] = []
#     lines.append("# Generated from Systre stdout (relaxed) -> CRYSTAL2")
#     lines.append(f"# Component: {comp.component_index}")
#     lines.append(f"# Systre-dimension: {comp.dimension}")
#     lines.append(f"# TD10: {comp.td10 or 'N/A'}")
#     lines.append("CRYSTAL")
#     lines.append(f"  NAME {comp.rcsr_name or name}")
#     lines.append(f"  GROUP {comp.ideal_space_group or 'P1'}")
#     lines.append(f"  CELL {a:.6f} {b:.6f} {c:.6f} {alpha:.6f} {beta:.6f} {gamma:.6f}")

#     for i, (x, y, z) in enumerate(node_coords, start=1):
#         lines.append(f"  NODE {i} {deg[i-1]}  {x:.6f} {y:.6f} {z:.6f}")

#     for (x1, y1, z1), (x2, y2, z2) in clean_edges:
#         lines.append(f"  EDGE  {x1:.6f} {y1:.6f} {z1:.6f}   {x2:.6f} {y2:.6f} {z2:.6f}")

#     for (x, y, z) in edge_centers:
#         lines.append(f"  EDGE_CENTER {x:.6f} {y:.6f} {z:.6f}")

#     lines.append("END")
#     out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

def write_crystal2_from_relaxed(out_path: Path, *, name: str, comp: RelaxedComponent) -> None:
    """
    Write a CRYSTAL2 CGD from a parsed relaxed component.

    Strategy (robust/general):
      - NODE lines come ONLY from Systre "Relaxed positions" (unique nodes).
      - NODE degree is taken from Systre "Coordination sequences" (first number), if available.
      - EDGE lines are written exactly as Systre printed them under "Edges:".
      - EDGE_CENTER lines are metadata (non-standard CGD token, but fine if your reader tolerates it).
    """

    # Sort nodes by label "1,2,3" if numeric; else lexicographic
    node_items = sorted(comp.nodes.items(), key=lambda kv: _sort_key(kv[0]))

    # Build degree map from coordination sequences (preferred)
    coord_deg: Dict[str, int] = getattr(comp, "coordination_number", {}) or {}

    # Fallback degree estimate from printed edges:
    # Count incidences where an edge endpoint equals a node coordinate (exact match).
    # This fallback is *best-effort* only; for many nets Systre prints translational endpoints
    # that won't match node coordinates exactly.
    fallback_deg_by_coord: Dict[Tuple[float, float, float], int] = {}
    for p, q in comp.edges:
        fallback_deg_by_coord[p] = fallback_deg_by_coord.get(p, 0) + 1
        fallback_deg_by_coord[q] = fallback_deg_by_coord.get(q, 0) + 1

    a, b, c, alpha, beta, gamma = comp.cell

    lines: List[str] = []
    lines.append("# Generated from Systre stdout (relaxed) -> CRYSTAL2")
    lines.append(f"# Component: {comp.component_index}")
    lines.append(f"# Systre-dimension: {comp.dimension}")
    lines.append(f"# TD10: {comp.td10 if comp.td10 is not None else 'N/A'}")
    if comp.rcsr_name:
        lines.append(f"# RCSR: {comp.rcsr_name}")

    lines.append("CRYSTAL")
    lines.append(f"  NAME {comp.rcsr_name or name}")
    lines.append(f"  GROUP {comp.ideal_space_group or 'P1'}")
    lines.append(f"  CELL {a:.6f} {b:.6f} {c:.6f} {alpha:.6f} {beta:.6f} {gamma:.6f}")

    # ---- NODES (general) ----
    for i, (label, (x, y, z)) in enumerate(node_items, start=1):
        deg_i = coord_deg.get(label)

        if deg_i is None:
            # best-effort fallback: look up the node coordinate in edge-incidence counts
            deg_i = fallback_deg_by_coord.get((x, y, z), 0)

        lines.append(f"  NODE {i} {deg_i}  {x:.6f} {y:.6f} {z:.6f}")

    # ---- EDGES (as printed by Systre) ----
    # Keep your de-dup if you want, but be careful: p<=q comparison is lexicographic and fine for tuples.
    seen = set()
    clean_edges: List[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = []
    for p, q in comp.edges:
        key = (p, q) if p <= q else (q, p)
        if key in seen:
            continue
        seen.add(key)
        clean_edges.append((p, q))

    for (x1, y1, z1), (x2, y2, z2) in clean_edges:
        lines.append(f"  EDGE  {x1:.6f} {y1:.6f} {z1:.6f}   {x2:.6f} {y2:.6f} {z2:.6f}")

    # ---- EDGE CENTERS (non-standard token; comment out if strict parser) ----
    for (x, y, z) in comp.edge_centers:
        lines.append(f"  # EDGE_CENTER  {x:.6f} {y:.6f} {z:.6f}")

    lines.append("END")
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
# =========================
# Fallback output writer
# =========================

def write_crystal2_fallback_from_periodic_graph(out_path: Path, *, name: str, periodic_graph_cgd: str) -> None:
    """
    Fallback: if Systre output doesn't include "relaxed" data, still produce a CRYSTAL2 CGD.

    Strategy: Convert the PERIODIC_GRAPH CGD text to a CRYSTAL2 wrapper.
    (This is not a "relaxed" cell, but it gives you a valid CRYSTAL2 CGD to keep pipeline working.)
    """
    lines_in = periodic_graph_cgd.splitlines()

    # Replace PERIODIC_GRAPH header with CRYSTAL
    out: List[str] = []
    out.append("# Fallback: CRYSTAL2 written from input PERIODIC_GRAPH (no relaxed output detected)")
    out.append("CRYSTAL")
    out.append(f"  NAME {name}")
    out.append("  GROUP P1")

    # Copy CELL/NODE/EDGE lines if present; otherwise just dump the payload
    for ln in lines_in:
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        if s.upper().startswith("PERIODIC_GRAPH"):
            continue
        if s.upper().startswith("NAME "):
            continue
        if s.upper().startswith("GROUP "):
            continue

        # Most CGD writers already output CELL/NODE/EDGE; keep them
        if s.upper().startswith(("CELL", "NODE", "EDGE")):
            out.append("  " + s)
        else:
            # keep other tokens (e.g. lattice vectors formats) indented as a best-effort
            out.append("  " + s)

    out.append("")
    out.append("END")
    out_path.write_text("\n".join(out) + "\n", encoding="utf-8")


# =========================
# CLI
# =========================

def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Build PERIODIC_GRAPH via mofstructure.generate_cgd, run Systre, save stdout, and write CRYSTAL2."
    )
    p.add_argument("input_structure", help="Input structure (e.g. CIF).")
    p.add_argument("--name", default="net")

    # Defaults are package-relative and will resolve from the installed package
    p.add_argument("--jar", default="bin/Systre-19.6.0.jar", help="Path to Systre JAR (or packaged resource).")
    p.add_argument("--arc", default="db/RCSRnets-2019-06-01.arc", help="Path to Systre ARC (or packaged resource).")

    p.add_argument("--java", default="java")
    p.add_argument("--xmx", default="1024m")
    p.add_argument("--timeout", type=int, default=60)

    p.add_argument("--systre-log", default="systre.cgd", help="Write Systre stdout here.")
    p.add_argument("--systre-err", default="systre.stderr", help="Write Systre stderr here.")
    p.add_argument("-o", "--output", default="systre_crystal2.cgd", help="Output CRYSTAL2 CGD.")
    p.add_argument("--keep-periodic-graph", default=None, help="Optional path to save the PERIODIC_GRAPH CGD fed to Systre.")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_argparser().parse_args(argv)

    inp = Path(args.input_structure)
    if not inp.exists():
        raise FileNotFoundError(inp)

    # 1) build PERIODIC_GRAPH with your  generate_cgd module
    topo = generate_cgd.TopologyExtractor(filename=str(inp))
    pg_text = topo.build_cgd(name=args.name)

    if args.keep_periodic_graph:
        Path(args.keep_periodic_graph).write_text(pg_text, encoding="utf-8")

    # write periodic graph to temp file for Systre
    fd, tmp_pg = tempfile.mkstemp(prefix="mofstructure_systre_", suffix=".cgd")
    os.close(fd)
    tmp_pg_path = Path(tmp_pg)
    tmp_pg_path.write_text(pg_text, encoding="utf-8")

    systre_log_path = Path(args.systre_log)
    systre_err_path = Path(args.systre_err)
    out_path = Path(args.output)

    try:
        # 2) run systre
        with systre_resources(args.jar, args.arc) as (jar, arc):
            print(f"Using Systre JAR: {jar}")
            print(f"Using Systre ARC: {arc}")

            proc = run_systre_on_cgd(
                tmp_pg_path,
                jar=jar,
                arc=arc,
                java=args.java,
                xmx=args.xmx,
                timeout_s=args.timeout,
            )

        # Always save outputs
        systre_log_path.write_text(proc.stdout or "", encoding="utf-8")
        systre_err_path.write_text(proc.stderr or "", encoding="utf-8")

        # If Systre failed, raise with stderr (but logs are saved)
        if proc.returncode != 0:
            raise RuntimeError(
                "Systre failed.\n"
                f"Return code: {proc.returncode}\n"
                f"Saved stdout to: {systre_log_path}\n"
                f"Saved stderr to: {systre_err_path}\n\n"
                f"Stderr:\n{proc.stderr}"
            )

        # 3) parse relaxed output (if present)
        stdout_text = proc.stdout or ""
        comps: List[RelaxedComponent] = []
        for idx, block in _split_component_blocks(stdout_text):
            comp = parse_relaxed_component(block, component_index=idx)
            if comp:
                comps.append(comp)

        if comps:
            comps.sort(key=lambda c: (len(c.edges), len(c.nodes)), reverse=True)
            best = comps[0]
            write_crystal2_from_relaxed(out_path, name=args.name, comp=best)
            print(f"Wrote CRYSTAL2 (from relaxed output): {out_path}")
        else:
            # 4) fallback: still produce a CRYSTAL2 CGD from the PERIODIC_GRAPH
            write_crystal2_fallback_from_periodic_graph(out_path, name=args.name, periodic_graph_cgd=pg_text)
            print(
                "WARNING: No 'relaxed' data detected in Systre stdout.\n"
                f" - Saved Systre stdout to: {systre_log_path}\n"
                f" - Saved Systre stderr to: {systre_err_path}\n"
                f" - Wrote fallback CRYSTAL2 from input PERIODIC_GRAPH: {out_path}\n"
                "If you expect relaxed output, open systre.cgd and check what Systre printed."
            )

        return 0

    finally:
        try:
            tmp_pg_path.unlink()
        except Exception:
            pass


if __name__ == "__main__":
    raise SystemExit(main())