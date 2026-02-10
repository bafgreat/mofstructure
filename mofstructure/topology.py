#!/usr/bin/env python3
"""
topology.py
===============

Centroid-based periodic topology extraction for porous crystals (MOFs/COFs) with CGD output.

This module provides a robust interface for generating **PERIODIC_GRAPH**
CGD files suitable for downstream topology identification (e.g., via Systre).

It wraps `mofstructure.mofdeconstructor` deconstruction routines and converts their
outputs (components + breaking pairs + optional region partitioning) into a periodic
graph representation.

-------------------------------------------------------------------------------
Conceptual overview
-------------------------------------------------------------------------------

Given a periodic atomic structure with a 3D unit cell, we define a *component graph*:

- **Nodes** represent molecular fragments / SBUs / ligands-clusters (depending on method).
  Each node corresponds to a *component*, i.e., a list of atom indices returned by
  `mofstructure.mofdeconstructor`.

- **Edges** are derived from *breaking pairs*:
  `breaking_pairs = [(i, j), ...]` are atom-index pairs across which connectivity was
  cut to separate components. Each (i, j) indicates that the component containing i
  is connected to the component containing j in the original structure.

Periodic shift (translation vector)
-----------------------------------

CGD edges require an integer lattice shift (sx, sy, sz) describing how the edge crosses
unit-cell boundaries. We compute this shift robustly using **fractional centroids**:

1. For each component, compute a periodic-consistent centroid in fractional coordinates.
   Because a component can straddle unit-cell boundaries, we unwrap each atom position
   relative to a reference atom within that component (minimum-image style in fractional
   space), then average, then wrap back into [0, 1).

2. For an edge from component u to v, let:
      delta = centroid_frac[v] - centroid_frac[u]
   The nearest-image convention suggests that the translation needed to bring v "closest"
   to u is:
      shift = -round(delta)
   where round is applied per component (np.rint). This yields integer shifts.

Edge multiplicity and deduplication
-----------------------------------

Breaking pairs can generate multiple edges between the same components. Some workflows
want to preserve multiplicity (coordination counts), while others (notably Systre) are
sensitive to duplicates. This module supports:

- no deduplication (preserve multiplicity)
- dedup exact parallel edges including shift (recommended for Systre)
- aggressive topological dedup ignoring shift (rarely recommended)

Region-target contraction (optional)
------------------------------------

`cgd_from_sbus_region_targets()` supports selecting the most connected *regions* and
contracting non-target components into effective edges between target components, while
reusing the already computed periodic edge shifts from the base component graph.

-------------------------------------------------------------------------------
Requirements
-------------------------------------------------------------------------------
- ASE (`ase`)
- mofstructure (`mofstructure.mofdeconstructor`)
- numpy

-------------------------------------------------------------------------------
Usage
-------------------------------------------------------------------------------

CLI:
    mofstructure_generate_gcd input.cif -o output.cgd --method all_node

Programmatic:
    from mofstructure.topology import TopologyExtractor

    topo = TopologyExtractor(filename="input.cif")
    cgd = topo.build_cgd(method="sbus", name="net")
    topo.write_cgd(cgd, "output.cgd")

Author: Dr. Dinga Wonanke
Status: production
"""

from __future__ import annotations

import argparse
import logging
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np
from ase.atoms import Atoms
from ase.io import read

from mofstructure import mofdeconstructor


logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


AtomIndex = Union[int, np.integer]
BreakingPair = Sequence[AtomIndex]
Component = Sequence[int]
Edge = Tuple[int, int, int, int, int]  # (u, v, sx, sy, sz)



def _component_centroid_periodic(frac_unwrapped: np.ndarray, comp: Sequence[int]) -> np.ndarray:
    """
    Compute a component centroid in fractional coordinates, robust to boundary wrapping.

    Parameters
    ----------
    frac_unwrapped
        Fractional coordinates from ASE with wrap=False, shape (n_atoms, 3).
    comp
        Atom indices belonging to the component.

    Returns
    -------
    centroid_frac
        Fractional centroid wrapped to [0, 1), shape (3,).
    """
    comp_i = [int(i) for i in comp]
    if not comp_i:
        raise ValueError("Empty component encountered while computing centroid.")

    f_ref = frac_unwrapped[comp_i[0]].copy()
    pts: List[np.ndarray] = []

    for idx in comp_i:
        fi = frac_unwrapped[idx].copy()
        # unwrap relative to reference in fractional space
        fi = fi - np.rint(fi - f_ref)
        pts.append(fi)

    c = np.mean(np.array(pts), axis=0)
    # wrap to [0, 1)
    c = c - np.floor(c)
    c[c >= 1.0] = 0.0
    return c


def cgd_from_components(
    atoms: Atoms,
    components: Sequence[Component],
    breaking_pairs: Sequence[BreakingPair],
    *,
    name: str = "net",
    spacegroup: str = "P1",
    dedup_parallel_edges: bool = True,
    dedup_mode: str = "shift",  # "none" | "shift" | "topological"
    sanity_tol: float = 1e-6,
    regions: Optional[Dict[int, List[int]]] = None,
) -> str:
    """
    Convert components + breaking pairs into CGD PERIODIC_GRAPH using centroid-based shifts.

    Parameters
    ----------
    atoms
        Periodic ASE Atoms object (must have 3D cell and pbc=True in all directions).
    components
        List of components. Each component is a list of atom indices.
    breaking_pairs
        List of (i, j) atom index pairs representing connections cut between components.
    name
        CGD graph ID.
    spacegroup
        Included for compatibility with external pipelines. (Not used in CGD content here.)
    dedup_parallel_edges
        Whether to deduplicate edges (recommended for Systre).
    dedup_mode
        - "none": keep all edges
        - "shift": remove exact duplicate (u,v,sx,sy,sz) edges (recommended)
        - "topological": remove duplicates ignoring shift (aggressive; usually avoid)
    sanity_tol
        Reserved for future geometric sanity checks (kept for API stability).
    regions
        Optional mapping: region_id -> list of component indices. Used for diagnostics only.

    Returns
    -------
    cgd_text
        CGD formatted string.

    Notes
    -----
    - Nodes are 1-based in CGD.
    - Shifts are integers describing translation in lattice coordinates.
    """
    _ = (spacegroup, sanity_tol)  # keep parameters for forward compatibility

    if atoms.cell is None or atoms.cell.rank < 3:
        raise ValueError("atoms must have a valid 3D unit cell.")
    if not np.all(atoms.pbc):
        raise ValueError("atoms.pbc must be True in all directions.")
    if not components:
        raise ValueError("components is empty.")

    n_components = len(components)

    # atom -> component
    atom_to_comp: Dict[int, int] = {}
    for ci, comp in enumerate(components):
        for a in comp:
            ai = int(a)
            if ai in atom_to_comp:
                raise ValueError(f"Atom index {ai} appears in multiple components.")
            atom_to_comp[ai] = ci

    # component centroids in fractional coordinates
    frac_unwrapped = atoms.get_scaled_positions(wrap=False)
    node_frac = np.zeros((n_components, 3), dtype=float)
    for ci, comp in enumerate(components):
        node_frac[ci] = _component_centroid_periodic(frac_unwrapped, comp)

    # region lookup (optional diagnostics)
    comp_to_region: Optional[Dict[int, int]] = None
    if regions is not None:
        comp_to_region = {}
        for rid, comp_ids in regions.items():
            for cid in comp_ids:
                comp_to_region[int(cid)] = int(rid)

    edges_raw: List[Edge] = []
    cross_region = same_region = unknown_region = 0

    for pair in breaking_pairs:
        if len(pair) != 2:
            raise ValueError(f"breaking_pairs element must be length-2, got: {pair}")

        i, j = int(pair[0]), int(pair[1])
        if i not in atom_to_comp or j not in atom_to_comp:
            raise ValueError(f"Breaking pair ({i},{j}) references atom not present in components.")

        u = atom_to_comp[i]
        v = atom_to_comp[j]
        if u == v:
            continue

        delta = node_frac[v] - node_frac[u]
        shift = -np.rint(delta).astype(int)
        sx, sy, sz = int(shift[0]), int(shift[1]), int(shift[2])

        if comp_to_region is not None:
            ru = comp_to_region.get(u, None)
            rv = comp_to_region.get(v, None)
            if ru is None or rv is None:
                unknown_region += 1
            elif ru == rv:
                same_region += 1
            else:
                cross_region += 1

        edges_raw.append((u, v, sx, sy, sz))

    # Deduplication
    edges_out: List[Edge] = edges_raw
    mode = dedup_mode.lower().strip()
    if dedup_parallel_edges:
        if mode not in {"none", "shift", "topological"}:
            raise ValueError("dedup_mode must be one of: 'none', 'shift', 'topological'")
        if mode != "none":
            seen = set()
            deduped: List[Edge] = []
            for (u, v, sx, sy, sz) in edges_out:
                if mode == "topological":
                    key = (min(u, v), max(u, v))
                else:
                    # canonical orientation in (u,v,shift)
                    if u < v:
                        key = (u, v, sx, sy, sz)
                    elif u > v:
                        key = (v, u, -sx, -sy, -sz)
                    else:
                        # self edge: canonicalize sign
                        key = (u, v, sx, sy, sz) if (sx, sy, sz) <= (0, 0, 0) else (u, v, -sx, -sy, -sz)

                if key in seen:
                    continue
                seen.add(key)
                deduped.append((u, v, sx, sy, sz))
            edges_out = deduped

    nonzero = sum(1 for (_, _, sx, sy, sz) in edges_out if (sx, sy, sz) != (0, 0, 0))

    def one_based(x: int) -> int:
        return x + 1

    # Render CGD
    lines: List[str] = []
    lines.append("# Generated by cgd_from_components (PERIODIC_GRAPH, centroid-based shifts)")
    lines.append(f"# Nodes (components): {n_components}")
    lines.append(f"# Breaking pairs: {len(breaking_pairs)}")
    lines.append(f"# Edges written: {len(edges_out)}")
    if comp_to_region is not None:
        lines.append(f"# Region connectivity: cross={cross_region}, same={same_region}, unknown={unknown_region}")
    lines.append(f"# Non-zero shifts: {nonzero}/{max(1, len(edges_out))}")
    lines.append(f"# Dedup: {dedup_parallel_edges} (mode={dedup_mode})")

    lines.append("PERIODIC_GRAPH")
    lines.append(f"ID {name}")
    lines.append("EDGES")
    for (u, v, sx, sy, sz) in edges_out:
        lines.append(f"  {one_based(u)} {one_based(v)} {sx} {sy} {sz}")
    lines.append("END")
    return "\n".join(lines) + "\n"


def cgd_from_sbus_region_targets(
    atoms: Atoms,
    components: Sequence[Component],
    breaking_pairs: Sequence[BreakingPair],
    regions: Dict[int, List[int]],
    *,
    name: str = "net",
    top_k_regions: Optional[int] = None,
    min_region_degree: Optional[int] = None,
    region_degree_percentile: Optional[float] = None,
    connect_mode: str = "clique",  # "clique" | "chain"
    dedup_edges_for_systre: bool = True,
) -> str:
    """
    Build PERIODIC_GRAPH where:
      - Nodes are components from the selected most-connected region(s)
      - All non-target components are contracted into effective edges between target nodes
      - Periodic shifts are inherited by reusing the base component graph from
        `cgd_from_components(dedup_mode="shift")`.

    Selection strategy (choose one):
      - top_k_regions: take K highest-scoring regions
      - min_region_degree: take all regions with score >= threshold
      - region_degree_percentile: take regions above a percentile
      - default: take all regions tied for the maximum score
    """
    if atoms.cell is None or atoms.cell.rank < 3:
        raise ValueError("atoms must have a valid 3D unit cell.")
    if not np.all(atoms.pbc):
        raise ValueError("atoms.pbc must be True in all directions.")
    if not regions:
        raise ValueError("regions is required for this method.")
    if not components:
        raise ValueError("components is empty.")

    n_components = len(components)

    # Build base periodic edges using the trusted routine, then parse EDGES back out.
    base_cgd = cgd_from_components(
        atoms=atoms,
        components=components,
        breaking_pairs=breaking_pairs,
        name="__base__",
        dedup_parallel_edges=True,
        dedup_mode="shift",
        regions=None,
    )

    base_edges: List[Edge] = []
    in_edges = False
    for line in base_cgd.splitlines():
        s = line.strip()
        if s == "EDGES":
            in_edges = True
            continue
        if s == "END":
            in_edges = False
        if in_edges and s and s[0].isdigit():
            parts = s.split()
            if len(parts) == 5:
                u1, v1, sx, sy, sz = map(int, parts)
                base_edges.append((u1 - 1, v1 - 1, sx, sy, sz))  # back to 0-based

    # adjacency + shift lookup
    adj: Dict[int, List[int]] = defaultdict(list)
    edge_shift: Dict[Tuple[int, int], List[Tuple[int, int, int]]] = defaultdict(list)
    for u, v, sx, sy, sz in base_edges:
        adj[u].append(v)
        adj[v].append(u)
        edge_shift[(u, v)].append((sx, sy, sz))
        edge_shift[(v, u)].append((-sx, -sy, -sz))

    # component -> region
    comp_to_region: Dict[int, int] = {}
    for rid, comp_ids in regions.items():
        for cid in comp_ids:
            comp_to_region[int(cid)] = int(rid)

    # score regions by incident base edges
    region_score = Counter()
    for u, v, *_ in base_edges:
        ru = comp_to_region.get(u, None)
        rv = comp_to_region.get(v, None)
        if ru is not None:
            region_score[ru] += 1
        if rv is not None:
            region_score[rv] += 1

    scores = {rid: region_score.get(rid, 0) for rid in regions.keys()}

    # choose target regions
    if top_k_regions is not None:
        target_regions = set([rid for rid, _ in sorted(scores.items(), key=lambda x: x[1], reverse=True)[:top_k_regions]])
    elif min_region_degree is not None:
        target_regions = set([rid for rid, s in scores.items() if s >= min_region_degree])
    elif region_degree_percentile is not None:
        vals = np.array(list(scores.values()), dtype=float)
        thr = float(np.percentile(vals, region_degree_percentile)) if len(vals) else 0.0
        target_regions = set([rid for rid, s in scores.items() if s >= thr])
    else:
        best = max(scores.values()) if scores else 0
        target_regions = set([rid for rid, s in scores.items() if s == best])

    if not target_regions:
        raise ValueError("No target regions selected. Relax thresholds or use top_k_regions.")

    target_components = sorted([c for c in range(n_components) if comp_to_region.get(c, None) in target_regions])
    if not target_components:
        raise ValueError("Target regions selected, but no components found in those regions.")

    node_id_of_comp = {c: i for i, c in enumerate(target_components)}
    edge_components = [c for c in range(n_components) if c not in node_id_of_comp]

    def best_shift(u: int, v: int) -> Optional[Tuple[int, int, int]]:
        lst = edge_shift.get((u, v), [])
        if not lst:
            return None
        return Counter(lst).most_common(1)[0][0]

    edges_out: List[Edge] = []
    mode = connect_mode.lower().strip()
    if mode not in {"clique", "chain"}:
        raise ValueError("connect_mode must be 'clique' or 'chain'")

    for ecomp in edge_components:
        touched = [nb for nb in adj.get(ecomp, []) if nb in node_id_of_comp]
        touched_u = sorted(set(touched))
        if len(touched_u) < 2:
            continue

        if mode == "chain":
            pairs = list(zip(touched_u[:-1], touched_u[1:]))
        else:
            pairs = [(touched_u[a], touched_u[b]) for a in range(len(touched_u)) for b in range(a + 1, len(touched_u))]

        for cu, cv in pairs:
            s1 = best_shift(cu, ecomp)
            s2 = best_shift(ecomp, cv)
            if s1 is None or s2 is None:
                continue

            sx = int(s1[0] + s2[0])
            sy = int(s1[1] + s2[1])
            sz = int(s1[2] + s2[2])

            u = node_id_of_comp[cu]
            v = node_id_of_comp[cv]
            if u == v and (sx, sy, sz) == (0, 0, 0):
                continue

            edges_out.append((u, v, sx, sy, sz))

    # dedup for Systre compatibility
    if dedup_edges_for_systre:
        seen = set()
        deduped: List[Edge] = []
        for (u, v, sx, sy, sz) in edges_out:
            if u < v:
                key = (u, v, sx, sy, sz)
            elif u > v:
                key = (v, u, -sx, -sy, -sz)
            else:
                key = (u, v, sx, sy, sz) if (sx, sy, sz) <= (0, 0, 0) else (u, v, -sx, -sy, -sz)
            if key in seen:
                continue
            seen.add(key)
            deduped.append((u, v, sx, sy, sz))
        edges_out = deduped

    nonzero = sum(1 for (_, _, sx, sy, sz) in edges_out if (sx, sy, sz) != (0, 0, 0))

    def one_based(x: int) -> int:
        return x + 1

    lines: List[str] = []
    lines.append("# Generated by cgd_from_sbus_region_targets (PERIODIC_GRAPH, contracted regions)")
    lines.append(f"# Target regions: {sorted(list(target_regions))}")
    lines.append("# Region scores: " + ", ".join([f"{rid}:{scores[rid]}" for rid in sorted(scores.keys())]))
    lines.append(f"# Nodes (target components): {len(target_components)} of {n_components} components")
    lines.append(f"# Edge-components contracted: {len(edge_components)}")
    lines.append(f"# Base edges: {len(base_edges)}")
    lines.append(f"# Contracted edges written: {len(edges_out)}")
    lines.append(f"# Non-zero shifts: {nonzero}/{max(1, len(edges_out))}")
    lines.append(f"# Dedup for Systre: {dedup_edges_for_systre}")

    lines.append("PERIODIC_GRAPH")
    lines.append(f"ID {name}")
    lines.append("EDGES")
    for (u, v, sx, sy, sz) in edges_out:
        lines.append(f"  {one_based(u)} {one_based(v)} {sx} {sy} {sz}")
    lines.append("END")
    return "\n".join(lines) + "\n"



@dataclass
class TopologyExtractor:
    """
    High-level façade for topology extraction and CGD writing.

    Parameters
    ----------
    ase_atoms
        ASE Atoms object (optional if `filename` is provided).
    filename
        Path to a structure file readable by ASE (e.g., CIF).
    """

    ase_atoms: Optional[Atoms] = None
    filename: Optional[str] = None

    def __post_init__(self) -> None:
        if self.ase_atoms is None:
            if not self.filename:
                raise ValueError("Provide either ase_atoms or filename.")
            self.ase_atoms = read(self.filename)

    @property
    def atoms(self) -> Atoms:
        assert self.ase_atoms is not None
        return self.ase_atoms

    def remove_guest(self) -> Atoms:
        """
        Remove unbound guests from a porous periodic structure.

        Returns
        -------
        guest_free_atoms
            ASE Atoms object with guest atoms removed.
        """
        idx = mofdeconstructor.remove_unbound_guest(self.atoms)
        return self.atoms[idx]

    def build_cgd(
        self,
        *,
        method: str = "sbus",
        name: str = "net",
        top_k_regions: int = 1,
        connect_mode: str = "clique",
        dedup_parallel_edges: bool = True,
        dedup_mode: str = "shift",
    ) -> str:
        """
        Build a CGD PERIODIC_GRAPH using one of the supported node-definition methods.

        Parameters
        ----------
        method
            - "sbus": nodes from `secondary_building_units`
            - "ligand_cluster": nodes from `ligands_and_metal_clusters`
            - "all_node": region-target contraction built on SBUs (your "all_node" mode)
        name
            CGD graph ID.
        top_k_regions
            Used for method="all_node": pick K most connected regions.
        connect_mode
            Used for method="all_node": "clique" or "chain".
        dedup_parallel_edges
            For "sbus" and "ligand_cluster": pass through to `cgd_from_components`.
        dedup_mode
            For "sbus" and "ligand_cluster": "none" | "shift" | "topological".

        Returns
        -------
        cgd_text
            CGD string.
        """
        guest_free = self.remove_guest()

        m = method.lower().strip()
        if m == "sbus":
            components, breaking_pairs, _porphyrin, regions = mofdeconstructor.secondary_building_units(guest_free)
            if not components:
                raise RuntimeError("Deconstruction failed: no components produced (secondary_building_units).")
            return cgd_from_components(
                atoms=guest_free,
                components=components,
                breaking_pairs=breaking_pairs,
                name=name,
                dedup_parallel_edges=dedup_parallel_edges,
                dedup_mode=dedup_mode,
                regions=regions,
            )

        if m in {"ligand_cluster", "ligands_cluster", "ligands_and_clusters"}:
            components, breaking_pairs, _porphyrin, regions = mofdeconstructor.ligands_and_metal_clusters(guest_free)
            if not components:
                raise RuntimeError("Deconstruction failed: no components produced (ligands_and_metal_clusters).")
            return cgd_from_components(
                atoms=guest_free,
                components=components,
                breaking_pairs=breaking_pairs,
                name=name,
                dedup_parallel_edges=dedup_parallel_edges,
                dedup_mode=dedup_mode,
                regions=regions,
            )

        if m in {"all_node", "region_targets", "sbus_region_targets"}:
            components, breaking_pairs, _porphyrin, regions = mofdeconstructor.secondary_building_units(guest_free)
            if not components:
                raise RuntimeError("Deconstruction failed: no components produced (secondary_building_units).")
            return cgd_from_sbus_region_targets(
                atoms=guest_free,
                components=components,
                breaking_pairs=breaking_pairs,
                regions=regions,
                name=name,
                top_k_regions=top_k_regions,
                connect_mode=connect_mode,
                dedup_edges_for_systre=True,  # critical for systre-like tools
            )

        raise ValueError("Unknown method. Choose from: 'sbus', 'all_node', 'ligand_cluster'.")

    @staticmethod
    def write_cgd(cgd_text: str, path: str) -> str:
        """Write CGD content to a file and return the output path."""
        with open(path, "w", encoding="utf-8") as f:
            f.write(cgd_text)
        return path



def build_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract periodic topology from a structure file and write CGD (PERIODIC_GRAPH)."
    )
    parser.add_argument("input_file", help="Path to input structure file (e.g., CIF).")
    parser.add_argument("-o", "--output", default="output.cgd", help="Path to output CGD file (default: output.cgd).")
    parser.add_argument(
        "--method",
        choices=["sbus", "all_node", "ligand_cluster"],
        default="sbus",
        help="Node-definition method.",
    )
    parser.add_argument("--name", default="net", help="CGD graph ID (default: net).")
    parser.add_argument(
        "--dedup",
        action="store_true",
        default=True,
        help="Deduplicate parallel edges (default: True).",
    )
    parser.add_argument(
        "--no-dedup",
        action="store_false",
        dest="dedup",
        help="Disable edge deduplication (keeps multiplicity).",
    )
    parser.add_argument(
        "--dedup-mode",
        choices=["none", "shift", "topological"],
        default="shift",
        help="Edge deduplication mode (default: shift).",
    )
    parser.add_argument(
        "--top-k-regions",
        type=int,
        default=1,
        help="For method=all_node: number of top regions to target (default: 1).",
    )
    parser.add_argument(
        "--connect-mode",
        choices=["clique", "chain"],
        default="clique",
        help="For method=all_node: how to connect nodes touched by a contracted component.",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_argparser()
    args = parser.parse_args(argv)

    topo = TopologyExtractor(filename=args.input_file)

    try:
        cgd_text = topo.build_cgd(
            method=args.method,
            name=args.name,
            top_k_regions=args.top_k_regions,
            connect_mode=args.connect_mode,
            dedup_parallel_edges=args.dedup,
            dedup_mode=args.dedup_mode,
        )
    except Exception as exc:
        logger.error("Failed to build CGD: %s", exc, exc_info=True)
        return 2

    TopologyExtractor.write_cgd(cgd_text, args.output)
    logger.info("CGD written to %s", args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
