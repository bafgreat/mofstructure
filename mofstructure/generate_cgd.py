#!/usr/bin/env python3
from __future__ import annotations
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"

import argparse
import logging
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
from ase.atoms import Atoms
from ase.io import read

from mofstructure import mofdeconstructor
from mofstructure.mofdeconstructor import transition_metals


logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )


AtomIndex = Union[int, np.integer]
BreakingPair = Sequence[Union[int, np.integer]]
Component = Sequence[int]
Edge = Tuple[int, int, int, int, int]  # (u, v, sx, sy, sz)


def remove_broken_bonds_from_breaking_pairs(
    atom_graph: Dict[int, List[int]],
    bond_offsets: Dict[Tuple[int, int], List[Tuple[int, int, int]]],
    breaking_pairs: Sequence[BreakingPair],
):
    '''
    Remove all bonds listed in `breaking_pairs` from an atom-level periodic graph.

    The first two entries of each breaking pair correspond to the atom indices
    of a bond that was cut during deconstruction. This function reconstructs
    the graph of bonds that remain after those cuts, while preserving the
    associated lattice offsets of the kept bonds.

    **parameters:**
        atom_graph: python dictionary
            Atom-level connectivity graph where each key is an atom index and the
            value is a list of neighbouring atom indices.

        bond_offsets: python dictionary
            Dictionary mapping (i, j) atom-index tuples to one or more periodic
            neighbour offsets, e.g. {(i, j): [(sx, sy, sz), ...]}.

        breaking_pairs: list
            List of broken bonds. Each entry is expected to be of the form
            [i, j] or [i, j, sx, sy, sz].

    **returns:**
        kept_graph: python dictionary
            Atom-level connectivity graph after removing the broken bonds.

        kept_offsets: python dictionary
            Dictionary of offsets corresponding to the bonds kept in the graph.
    '''
    broken = set()
    for pair in breaking_pairs:
        if len(pair) < 2:
            continue
        i, j = int(pair[0]), int(pair[1])
        broken.add((i, j))
        broken.add((j, i))

    kept_graph = {}
    kept_offsets = defaultdict(list)

    for i, neighs in atom_graph.items():
        keep = []
        for j in neighs:
            j = int(j)
            if (int(i), j) in broken:
                continue
            keep.append(j)
            for s in bond_offsets.get((int(i), j), []):
                kept_offsets[(int(i), j)].append(tuple(int(x) for x in s))
        kept_graph[int(i)] = keep

    return kept_graph, kept_offsets


def component_atom_images(
    components: Sequence[Component],
    kept_graph: Dict[int, List[int]],
    kept_offsets: Dict[Tuple[int, int], List[Tuple[int, int, int]]],
):
    '''
    Assign integer lattice image vectors to atoms inside each connected component.

    After the bonds listed in `breaking_pairs` are removed, the structure is split
    into final connected components. Some of these components may still span the
    periodic unit cell. This function unwraps each component by traversing the kept
    internal bonds and assigning an integer image vector (ix, iy, iz) to each atom.

    These image vectors are required to compute the true translation between
    components connected by a broken bond.

    **parameters:**
        components: list
            List of connected components, where each component is a list of atom indices.

        kept_graph: python dictionary
            Atom-level connectivity graph after removing broken bonds.

        kept_offsets: python dictionary
            Dictionary mapping kept bonds to their periodic offsets.

    **returns:**
        comp_images: python dictionary
            Dictionary of the form:
                comp_images[component_id][atom_index] = (ix, iy, iz)
    '''
    comp_images = {}

    for cid, comp in enumerate(components):
        comp_set = set(int(a) for a in comp)
        if not comp_set:
            comp_images[cid] = {}
            continue

        root = next(iter(comp_set))
        images = {root: (0, 0, 0)}
        stack = [root]

        while stack:
            a = stack.pop()
            ia = np.array(images[a], dtype=int)

            for b in kept_graph.get(a, []):
                b = int(b)
                if b not in comp_set:
                    continue

                shifts = kept_offsets.get((a, b), [])
                if not shifts:
                    continue

                # Use the first stored offset for this kept bond.
                s = np.array(shifts[0], dtype=int)
                ib = tuple((ia + s).tolist())

                if b not in images:
                    images[b] = ib
                    stack.append(b)

        for a in comp_set:
            images.setdefault(a, (0, 0, 0))

        comp_images[cid] = images

    return comp_images


def base_edges_with_shifts(
    atoms: Atoms,
    components: Sequence[Component],
    breaking_pairs: Sequence[BreakingPair],
):
    '''
    Build the base component graph with periodic shifts.

    A broken pair returned by `mofdeconstructor` contains the indices of two atoms
    that were bonded before deconstruction, and optionally the local bond offset
    between them. However, this local bond offset alone is not sufficient to
    determine the true translation between the final post-cut components.

    This function therefore:
        1) reconstructs the kept-bond graph,
        2) unwraps each final component using the kept periodic bonds,
        3) computes the translation between components using:

            T(u -> v) = img_u[i] - img_v[j] + s_ij

    where:
        - i belongs to component u
        - j belongs to component v
        - img_u[i] is the lattice image of atom i in component u
        - img_v[j] is the lattice image of atom j in component v
        - s_ij is the bond offset stored in breaking_pairs

    **parameters:**
        atoms: ASE atoms object

        components: list
            List of final connected components after bond breaking.

        breaking_pairs: list
            List of broken bonds returned by the deconstructor.
            Each entry is typically [i, j, sx, sy, sz].

    **returns:**
        edges_out: list
            List of component-level edges of the form:
                (u, v, sx, sy, sz)
            where u and v are component indices and (sx, sy, sz) is the
            corresponding lattice shift.
    '''
    if atoms.cell is None or atoms.cell.rank < 3:
        raise ValueError("atoms must have a valid 3D unit cell.")
    if not np.all(atoms.pbc):
        raise ValueError("atoms.pbc must be True in all directions.")
    if not components:
        raise ValueError("components is empty.")

    atom_to_comp = {}
    for ci, comp in enumerate(components):
        for a in comp:
            ai = int(a)
            if ai in atom_to_comp:
                raise ValueError(f"Atom index {ai} appears in multiple components.")
            atom_to_comp[ai] = ci

    atom_graph, _, bond_offsets = mofdeconstructor.compute_ase_neighbour_with_offsets(atoms)

    kept_graph, kept_offsets = remove_broken_bonds_from_breaking_pairs(
        atom_graph, bond_offsets, breaking_pairs
    )

    comp_images = component_atom_images(components, kept_graph, kept_offsets)

    edges_raw = []

    for pair in breaking_pairs:
        if len(pair) < 2:
            continue

        i = int(pair[0])
        j = int(pair[1])

        if i not in atom_to_comp or j not in atom_to_comp:
            continue

        u = atom_to_comp[i]
        v = atom_to_comp[j]

        if u == v:
            continue

        sij = np.array([
            int(pair[2]) if len(pair) > 2 else 0,
            int(pair[3]) if len(pair) > 3 else 0,
            int(pair[4]) if len(pair) > 4 else 0,
        ], dtype=int)

        iu = np.array(comp_images[u].get(i, (0, 0, 0)), dtype=int)
        iv = np.array(comp_images[v].get(j, (0, 0, 0)), dtype=int)

        tuv = iu - iv + sij
        sx, sy, sz = int(tuv[0]), int(tuv[1]), int(tuv[2])

        edges_raw.append((u, v, sx, sy, sz))

    seen = set()
    edges_out = []

    for (u, v, sx, sy, sz) in edges_raw:
        if u < v:
            key = (u, v, sx, sy, sz)
            out = (u, v, sx, sy, sz)
        elif u > v:
            key = (v, u, -sx, -sy, -sz)
            out = (v, u, -sx, -sy, -sz)
        else:
            if (sx, sy, sz) > (0, 0, 0):
                sx, sy, sz = -sx, -sy, -sz
            key = (u, v, sx, sy, sz)
            out = (u, v, sx, sy, sz)

        if key in seen:
            continue
        seen.add(key)
        edges_out.append(out)

    return edges_out


def component_has_transition_metal_not_in_porphyrin(
    atoms: Atoms,
    comp: Sequence[int],
    porphyrin_atoms: Optional[Sequence[int]],
    tm_symbols: Sequence[str],
):
    '''
    Check whether a connected component contains a transition metal that is not
    part of a porphyrin-like motif.

    **parameters:**
        atoms: ASE atoms object

        comp: list
            List of atom indices defining a connected component.

        porphyrin_atoms: list or None
            List of atom indices identified as belonging to a porphyrin motif.

        tm_symbols: list
            List of transition-metal symbols.

    **returns:**
        bool
    '''
    porph_set = set(int(i) for i in (porphyrin_atoms or []))
    tm_set = set(str(s) for s in tm_symbols)

    for ai in comp:
        a = int(ai)
        if atoms[a].symbol in tm_set and a not in porph_set:
            return True
    return False


def regions_with_metal(
    atoms: Atoms,
    components: Sequence[Component],
    regions: Dict[int, List[int]],
    porphyrin_atoms: Optional[Sequence[int]],
):
    '''
    Identify regions that contain at least one metal-bearing component.

    In this module, a region corresponds to a group of equivalent components
    returned by the deconstructor. This function selects the regions that contain
    at least one connected component with a transition metal atom not assigned
    to a porphyrin motif.

    **parameters:**
        atoms: ASE atoms object

        components: list
            List of connected components.

        regions: python dictionary
            Mapping of region_id -> list of component indices.

        porphyrin_atoms: list or None
            Atom indices associated with porphyrin-like motifs.

    **returns:**
        list
            Sorted list of region IDs containing metal-bearing components.
    '''
    tm = transition_metals()
    target = set()

    for rid, comp_ids in regions.items():
        for cid in comp_ids:
            c = int(cid)
            if component_has_transition_metal_not_in_porphyrin(
                atoms, components[c], porphyrin_atoms, tm
            ):
                target.add(int(rid))
                break

    return sorted(target)


def cgd_from_region_targets(
    atoms: Atoms,
    components: Sequence[Component],
    breaking_pairs: Sequence[BreakingPair],
    regions: Dict[int, List[int]],
    *,
    target_regions: Sequence[int],
    name: str = "net",
    connect_mode: str = "clique",
    dedup_edges_for_systre: bool = True,
):
    '''
    Build a CGD PERIODIC_GRAPH by contracting non-target components into edges.

    The logic is:
        1) Build the base component graph with periodic shifts.
        2) Select all components belonging to `target_regions` as CGD nodes.
        3) Treat all remaining components as edge-components.
        4) Contract each edge-component into effective edges between node-components.

    The contraction is incidence-based, meaning that different periodic incidences
    between a linker-like component and node-like components are preserved. This is
    essential for recovering the correct 3D net of pillared or otherwise periodic
    frameworks.

    **parameters:**
        atoms: ASE atoms object

        components: list
            List of connected components from the deconstructor.

        breaking_pairs: list
            List of broken atom pairs from the deconstructor.

        regions: python dictionary
            Mapping of region_id -> list of component indices.

        target_regions: list
            Region IDs that should be treated as node regions.

        name: str
            CGD graph ID.

        connect_mode: str
            Either "clique" or "chain".
            - "clique": all incidences are mutually connected
            - "chain": incidences are connected in sorted order

        dedup_edges_for_systre: bool
            If True, remove exact duplicate periodic edges.

    **returns:**
        cgd_text: str
            CGD PERIODIC_GRAPH content as a string.
    '''
    if not regions:
        raise ValueError("regions is required for contraction methods.")
    if not target_regions:
        raise ValueError("target_regions is empty.")
    if not components:
        raise ValueError("components is empty.")

    n_components = len(components)
    target_regions = set(int(r) for r in target_regions)

    base_edges = base_edges_with_shifts(atoms, components, breaking_pairs)

    adj = defaultdict(list)
    edge_shift = defaultdict(list)

    for u, v, sx, sy, sz in base_edges:
        adj[u].append(v)
        adj[v].append(u)
        edge_shift[(u, v)].append((sx, sy, sz))
        edge_shift[(v, u)].append((-sx, -sy, -sz))

    comp_to_region = {}
    for rid, comp_ids in regions.items():
        for cid in comp_ids:
            comp_to_region[int(cid)] = int(rid)

    target_components = sorted(
        [c for c in range(n_components) if comp_to_region.get(c, None) in target_regions]
    )
    if not target_components:
        raise ValueError(f"Target regions {sorted(target_regions)} contain no components.")

    node_id_of_comp = {c: i for i, c in enumerate(target_components)}
    edge_components = [c for c in range(n_components) if c not in node_id_of_comp]

    mode = connect_mode.lower().strip()
    if mode not in {"clique", "chain"}:
        raise ValueError("connect_mode must be 'clique' or 'chain'")

    edges_out = []

    for ecomp in edge_components:
        incidences = []

        for nb in adj.get(ecomp, []):
            if nb not in node_id_of_comp:
                continue
            for s in edge_shift.get((nb, ecomp), []):
                incidences.append((nb, (int(s[0]), int(s[1]), int(s[2]))))

        if len(incidences) < 2:
            continue

        if mode == "chain":
            incidences = sorted(
                incidences,
                key=lambda x: (x[0], x[1][0], x[1][1], x[1][2])
            )
            incidence_pairs = list(zip(incidences[:-1], incidences[1:]))
        else:
            incidence_pairs = [
                (incidences[a], incidences[b])
                for a in range(len(incidences))
                for b in range(a + 1, len(incidences))
            ]

        for (cu, su), (cv, sv) in incidence_pairs:
            u = node_id_of_comp[cu]
            v = node_id_of_comp[cv]

            sx = int(su[0] - sv[0])
            sy = int(su[1] - sv[1])
            sz = int(su[2] - sv[2])

            if u == v and (sx, sy, sz) == (0, 0, 0):
                continue

            edges_out.append((u, v, sx, sy, sz))

    if dedup_edges_for_systre:
        seen = set()
        deduped = []

        for (u, v, sx, sy, sz) in edges_out:
            if u < v:
                key = (u, v, sx, sy, sz)
                out = (u, v, sx, sy, sz)
            elif u > v:
                key = (v, u, -sx, -sy, -sz)
                out = (v, u, -sx, -sy, -sz)
            else:
                if (sx, sy, sz) <= (0, 0, 0):
                    key = (u, v, sx, sy, sz)
                    out = (u, v, sx, sy, sz)
                else:
                    key = (u, v, -sx, -sy, -sz)
                    out = (u, v, -sx, -sy, -sz)

            if key in seen:
                continue
            seen.add(key)
            deduped.append(out)

        edges_out = deduped

    nonzero = sum(1 for (_, _, sx, sy, sz) in edges_out if (sx, sy, sz) != (0, 0, 0))

    def one_based(x):
        return x + 1

    lines = []
    lines.append("# Generated by mofstructure.topology (PERIODIC_GRAPH, region contraction)")
    lines.append(f"# Target regions: {sorted(list(target_regions))}")
    lines.append(f"# Nodes (target components): {len(target_components)}")
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
    '''
    A high-level class for extracting a periodic topology from a porous crystal.

    The class accepts either an ASE atoms object or a structure filename readable
    by ASE. It removes unbound guests, performs deconstruction using one of the
    supported methods, selects node regions, contracts linker-like regions, and
    writes the resulting CGD PERIODIC_GRAPH.

    **parameters:**
        ase_atoms: ASE atoms object, optional

        filename: str, optional
            Path to an input structure file readable by ASE.

    **methods:**
        remove_guest():
            Remove unbound guests from the framework.

        build_cgd():
            Build a CGD periodic graph using one of the supported topology modes.

        write_cgd():
            Write CGD content to file.
    '''
    ase_atoms: Optional[Atoms] = None
    filename: Optional[str] = None

    def __post_init__(self):
        if self.ase_atoms is None:
            if not self.filename:
                raise ValueError("Provide either ase_atoms or filename.")
            self.ase_atoms = read(self.filename)

    @property
    def atoms(self):
        return self.ase_atoms

    def remove_guest(self):
        '''
        Remove unbound guests from a porous periodic structure.

        **returns:**
            ASE atoms object
                Guest-free framework.
        '''
        idx = mofdeconstructor.remove_unbound_guest(self.atoms)
        return self.atoms[idx]

    def build_cgd(
        self,
        *,
        method: str = "sbus",
        name: str = "net",
        top_k_regions: int = 1,
        connect_mode: str = "clique",
    ):
        '''
        Build a CGD PERIODIC_GRAPH using one of the supported methods.

        **parameters:**
            method: str
                One of:
                    - "sbus"
                    - "ligand_cluster"
                    - "all_node"

            name: str
                CGD graph ID.

            top_k_regions: int
                Only used for method="all_node". Number of most-connected regions
                to keep as node regions.

            connect_mode: str
                Either "clique" or "chain".

        **returns:**
            cgd_text: str
                CGD PERIODIC_GRAPH content.
        '''
        guest_free = self.remove_guest()
        method = method.lower().strip()

        logger.info(
            "Building CGD with method='%s', name='%s', top_k_regions=%s, connect_mode='%s'",
            method, name, top_k_regions, connect_mode
        )

        if method == "sbus":
            components, _, porphyrin, regions, breaking_pairs = mofdeconstructor.secondary_building_units(guest_free)

            if not components:
                raise RuntimeError("Deconstruction failed: no components (secondary_building_units).")
            if not regions:
                raise RuntimeError("Method 'sbus' requires regions, but regions is empty.")

            target_regions = regions_with_metal(guest_free, components, regions, porphyrin)
            if not target_regions:
                raise RuntimeError("No metal-containing regions found (excluding porphyrin metals) for method='sbus'.")

            return cgd_from_region_targets(
                atoms=guest_free,
                components=components,
                breaking_pairs=breaking_pairs,
                regions=regions,
                target_regions=target_regions,
                name=name,
                connect_mode=connect_mode,
                dedup_edges_for_systre=True,
            )

        elif method == "ligand_cluster":
            components, _, porphyrin, regions, breaking_pairs = mofdeconstructor.ligands_and_metal_clusters(guest_free)

            if not components:
                raise RuntimeError("Deconstruction failed: no components (ligands_and_metal_clusters).")
            if not regions:
                raise RuntimeError("Method 'ligand_cluster' requires regions, but regions is empty.")

            target_regions = regions_with_metal(guest_free, components, regions, porphyrin)
            if not target_regions:
                raise RuntimeError("No metal-containing regions found (excluding porphyrin metals) for method='ligand_cluster'.")

            return cgd_from_region_targets(
                atoms=guest_free,
                components=components,
                breaking_pairs=breaking_pairs,
                regions=regions,
                target_regions=target_regions,
                name=name,
                connect_mode=connect_mode,
                dedup_edges_for_systre=True,
            )

        elif method == "all_node":
            components, _, _, regions, breaking_pairs = mofdeconstructor.secondary_building_units(guest_free)

            if not components:
                raise RuntimeError("Deconstruction failed: no components (secondary_building_units).")
            if not regions:
                raise RuntimeError("Method 'all_node' requires regions, but regions is empty.")

            base_edges = base_edges_with_shifts(guest_free, components, breaking_pairs)

            comp_to_region = {}
            for rid, comp_ids in regions.items():
                for cid in comp_ids:
                    comp_to_region[int(cid)] = int(rid)

            region_score = Counter()
            for u, v, *_ in base_edges:
                ru = comp_to_region.get(u, None)
                rv = comp_to_region.get(v, None)
                if ru is not None:
                    region_score[ru] += 1
                if rv is not None:
                    region_score[rv] += 1

            ranked = [rid for rid, _ in sorted(region_score.items(), key=lambda x: x[1], reverse=True)]
            target_regions = ranked[:max(1, int(top_k_regions))]

            if not target_regions:
                raise RuntimeError("Method 'all_node' could not select any target regions.")

            return cgd_from_region_targets(
                atoms=guest_free,
                components=components,
                breaking_pairs=breaking_pairs,
                regions=regions,
                target_regions=target_regions,
                name=name,
                connect_mode=connect_mode,
                dedup_edges_for_systre=True,
            )

        else:
            raise ValueError("Unknown method. Choose from: 'sbus', 'ligand_cluster', 'all_node'.")

    @staticmethod
    def write_cgd(cgd_text, path):
        '''
        Write CGD content to file.

        **parameters:**
            cgd_text: str
                CGD content.

            path: str
                Output file path.

        **returns:**
            str
                Output path.
        '''
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(cgd_text)
        return path


def build_argparser():
    '''
    Build the command-line argument parser.

    **returns:**
        argparse.ArgumentParser
    '''
    parser = argparse.ArgumentParser(
        description="Extract periodic topology from a structure file and write CGD (PERIODIC_GRAPH)."
    )
    parser.add_argument("input_file", help="Path to input structure file (e.g. CIF).")
    parser.add_argument(
        "-o", "--output",
        default="output.cgd",
        help="Path to output CGD file (default: output.cgd)."
    )
    parser.add_argument(
        "--method",
        choices=["sbus", "ligand_cluster", "all_node"],
        default="sbus",
        help="Topology extraction method."
    )
    parser.add_argument(
        "--name",
        default="net",
        help="CGD graph ID (default: net)."
    )
    parser.add_argument(
        "--top-k-regions",
        type=int,
        default=1,
        help="For method='all_node': number of top connected regions to keep as nodes."
    )
    parser.add_argument(
        "--connect-mode",
        choices=["clique", "chain"],
        default="clique",
        help="Contraction connectivity mode."
    )
    return parser


def main(argv: Optional[Sequence[str]] = None):
    '''
    Command-line entry point.

    **parameters:**
        argv: list, optional

    **returns:**
        int
            Exit code.
    '''
    parser = build_argparser()
    args = parser.parse_args(argv)

    topo = TopologyExtractor(filename=args.input_file)

    try:
        cgd_text = topo.build_cgd(
            method=args.method,
            name=args.name,
            top_k_regions=args.top_k_regions,
            connect_mode=args.connect_mode,
        )
    except Exception as exc:
        logger.error("Failed to build CGD: %s", exc, exc_info=True)
        return 2

    TopologyExtractor.write_cgd(cgd_text, args.output)
    logger.info("CGD written to %s", args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())