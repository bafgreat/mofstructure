#!/usr/bin/env python3
'''
Construction of CGD nets from crystal structures.

Systre works on nets rather than atoms, so a framework has to be reduced to
vertices and edges first. The node definition decides what the net describes,
and three are available: sbus places a vertex at each secondary building unit,
all_node keeps every branch point (splitting rod SBUs into their atoms), and
single_node coarsens all_node by merging each organic group to one vertex. The
last two reproduce CrystalNets' AllNodes and SingleNodes. The same framework can
give different nets under each, which is expected rather than an error.
'''
from __future__ import annotations
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"

import argparse
import logging
from collections import defaultdict
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


def _canonical_translation(vec: Tuple[int, int, int]):
    '''
    Reduce a lattice translation to a primitive, sign-canonical form so that a
    vector and its negative, and any integer multiple, all map to one key.

    A rod running along c is discovered as (0, 0, 1); the same rod found through
    a longer path may surface as (0, 0, 2). Dividing by the gcd collapses both
    to (0, 0, 1), and forcing the first non-zero component positive collapses
    (0, 0, 1) and (0, 0, -1) together.
    '''
    x, y, z = (int(vec[0]), int(vec[1]), int(vec[2]))
    if (x, y, z) == (0, 0, 0):
        return None

    g = np.gcd.reduce([abs(x), abs(y), abs(z)])
    if g > 1:
        x, y, z = x // g, y // g, z // g

    for component in (x, y, z):
        if component != 0:
            if component < 0:
                x, y, z = -x, -y, -z
            break
    return (x, y, z)


def component_self_translations(
    components: Sequence[Component],
    kept_graph: Dict[int, List[int]],
    kept_offsets: Dict[Tuple[int, int], List[Tuple[int, int, int]]],
):
    '''
    Find the lattice translations along which each component is periodic within
    itself.

    A discrete building unit (a paddlewheel, a Zr6 cluster) is finite: unwrapping
    its internal bonds assigns every atom one consistent image, so it has no
    self-translation. A rod-shaped SBU is an infinite chain, so following its
    internal bonds eventually returns to an atom already seen but in a different
    periodic image. That image difference is a translation of the rod, and it is
    exactly the connectivity that keeps a rod framework three-dimensional.

    ``component_atom_images`` discards these back-edges; this recovers them. Each
    translation is reduced to a primitive, sign-canonical vector, so a rod yields
    one vector and a discrete cluster yields none.

    **parameters:**
        components: list
            Connected components as lists of atom indices.

        kept_graph: python dictionary
            Atom-level graph after the broken bonds were removed.

        kept_offsets: python dictionary
            Periodic offsets of the kept bonds.

    **returns:**
        python dictionary
            Mapping of component id -> sorted list of primitive translations
            (sx, sy, sz). Empty for finite components.
    '''
    self_translations: Dict[int, List[Tuple[int, int, int]]] = {}

    for cid, comp in enumerate(components):
        comp_set = set(int(a) for a in comp)
        if not comp_set:
            self_translations[cid] = []
            continue

        root = next(iter(comp_set))
        images = {root: (0, 0, 0)}
        stack = [root]
        found = set()

        while stack:
            a = stack.pop()
            ia = np.array(images[a], dtype=int)

            for b in kept_graph.get(a, []):
                b = int(b)
                if b not in comp_set:
                    continue

                for s in kept_offsets.get((a, b), []):
                    ib = tuple((ia + np.array(s, dtype=int)).tolist())
                    if b not in images:
                        images[b] = ib
                        stack.append(b)
                    else:
                        drift = tuple(
                            int(ib[k] - images[b][k]) for k in range(3)
                        )
                        reduced = _canonical_translation(drift)
                        if reduced is not None:
                            found.add(reduced)

        self_translations[cid] = sorted(found)

    return self_translations


def kept_bond_graph(atoms, breaking_pairs):
    '''
    Build the atom-level graph that remains after the deconstruction cuts.

    Wraps the neighbour search and bond removal so a caller that needs the kept
    graph for more than one purpose (edges between components and each
    component's own periodicity, say) can compute it once and pass it around.

    **parameters:**
        atoms: ASE atoms object
        breaking_pairs: broken bonds from the deconstructor

    **returns:**
        (kept_graph, kept_offsets)
    '''
    atom_graph, _, bond_offsets = \
        mofdeconstructor.compute_ase_neighbour_with_offsets(atoms)
    return remove_broken_bonds_from_breaking_pairs(
        atom_graph, bond_offsets, breaking_pairs
    )


def base_edges_with_shifts(
    atoms: Atoms,
    components: Sequence[Component],
    breaking_pairs: Sequence[BreakingPair],
    kept_graph: Optional[Dict[int, List[int]]] = None,
    kept_offsets: Optional[Dict[Tuple[int, int], List[Tuple[int, int, int]]]] = None,
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
            List of component-level edges of the form::

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

    if kept_graph is None or kept_offsets is None:
        kept_graph, kept_offsets = kept_bond_graph(atoms, breaking_pairs)

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

    # The kept-bond graph is needed both for the edges between components and,
    # below, for each component's own periodicity, so compute it once and share.
    kept_graph, kept_offsets = kept_bond_graph(atoms, breaking_pairs)
    base_edges = base_edges_with_shifts(
        atoms, components, breaking_pairs, kept_graph, kept_offsets
    )

    # Rod SBUs are periodic within themselves, and that periodicity is what keeps
    # a rod framework three-dimensional. base_edges only carries connections
    # between components, so the intra-rod translations are recovered separately
    # and emitted as self-edges below. Discrete clusters return nothing here.
    self_translations = component_self_translations(
        components, kept_graph, kept_offsets
    )

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

    edges_out = []

    # linker nodes are numbered after the metal (target) nodes
    next_linker_node = len(target_components)

    for ecomp in edge_components:
        incidences = []

        for nb in adj.get(ecomp, []):
            if nb not in node_id_of_comp:
                continue
            for s in edge_shift.get((nb, ecomp), []):
                node = node_id_of_comp[nb]
                incidences.append((node, (int(s[0]), int(s[1]), int(s[2]))))

        # A linker binding the same metal node through more than one atom (a
        # chelating carboxylate, for instance) records the incidence twice with
        # the same shift. Collapse those so the linker's connectivity reflects
        # how many distinct metal sites it actually bridges.
        incidences = sorted(set(incidences))

        if len(incidences) < 2:
            continue

        if len(incidences) == 2:
            # A ditopic linker is a plain connection between two metal nodes, so
            # it contracts to a single edge. Keeping it as a 2-connected node
            # would only subdivide that edge, giving the same topology.
            (u, su), (v, sv) = incidences
            sx = int(su[0] - sv[0])
            sy = int(su[1] - sv[1])
            sz = int(su[2] - sv[2])
            if not (u == v and (sx, sy, sz) == (0, 0, 0)):
                edges_out.append((u, v, sx, sy, sz))
            continue

        # A polytopic linker (tritopic BTC and higher) is a branch point, so it
        # must be its own node joined to each metal it bridges. Contracting it
        # into a clique instead would inflate every metal's connectivity and
        # give the wrong net -- this is what turned HKUST-1 into reo (8-c) when
        # it should be tbo (a 3,4-connected net).
        lnode = next_linker_node
        next_linker_node += 1
        for (u, su) in incidences:
            edges_out.append((u, lnode, int(su[0]), int(su[1]), int(su[2])))

    # give each rod-shaped node its intra-chain connectivity back
    for comp_id, node_id in node_id_of_comp.items():
        for (sx, sy, sz) in self_translations.get(comp_id, []):
            edges_out.append((node_id, node_id, sx, sy, sz))

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


def dedup_periodic_edges(edges):
    '''
    Remove duplicate undirected periodic edges. An edge (u, v, s) and its reverse
    (v, u, -s) are the same edge; a self-edge (u, u, s) and (u, u, -s) likewise.
    '''
    seen = set()
    out = []
    for u, v, sx, sy, sz in edges:
        neg = (-sx, -sy, -sz)
        if (u, (sx, sy, sz)) <= (v, neg):
            key = (u, v, sx, sy, sz)
        else:
            key = (v, u, neg[0], neg[1], neg[2])
        if key[0] == key[1] and key[2:] == (0, 0, 0):
            continue
        if key in seen:
            continue
        seen.add(key)
        out.append(key)
    return out


def periodic_graph_cgd(edges, name):
    '''
    Format a list of (u, v, sx, sy, sz) edges (0-based nodes) as a CGD
    PERIODIC_GRAPH that Systre can read.
    '''
    lines = [
        "# Generated by mofstructure.topology (PERIODIC_GRAPH, all-node)",
        "PERIODIC_GRAPH",
        f"ID {name}",
        "EDGES",
    ]
    for u, v, sx, sy, sz in edges:
        lines.append(f"  {u + 1} {v + 1} {sx} {sy} {sz}")
    lines.append("END")
    return "\n".join(lines) + "\n"


def cgd_all_nodes(
    atoms: Atoms,
    components: Sequence[Component],
    breaking_pairs: Sequence[BreakingPair],
    regions: Dict[int, List[int]],
    *,
    target_regions: Sequence[int],
    name: str = "net",
):
    '''
    Build the all-node CGD, where a rod SBU keeps its atoms as separate nodes.

    A discrete metal cluster collapses to one node, exactly as
    ``cgd_from_region_targets`` does, so the net for a discrete-SBU framework is
    unchanged. A rod SBU is an infinite chain, so collapsing it to one node
    throws away the chain and drops the framework to a lower-dimensional net. It
    is instead split: each metal atom and each bridging carboxyl carbon becomes
    a node, and the oxygen atoms between them contract to edges. This recovers
    the true net -- MIL-53 gives rna rather than the collapsed pcu.

    **parameters:**
        atoms: ASE atoms object

        components: connected components from the deconstructor

        breaking_pairs: broken bonds from the deconstructor

        regions: mapping region_id -> component ids

        target_regions: region ids treated as node (metal) regions

        name: CGD graph ID

    **returns:**
        cgd_text: str
    '''
    edges, _ = _all_node_graph(
        atoms, components, breaking_pairs, regions, target_regions
    )
    return periodic_graph_cgd(dedup_periodic_edges(edges), name)


def cgd_single_nodes(
    atoms: Atoms,
    components: Sequence[Component],
    breaking_pairs: Sequence[BreakingPair],
    regions: Dict[int, List[int]],
    *,
    target_regions: Sequence[int],
    name: str = "net",
):
    '''
    Build the single-node CGD, a coarsening of the all-node net.

    Starting from the all-node net, every connected group of organic vertices
    (the carboxyl-carbon and linker nodes) is merged into a single vertex, while
    the metal vertices are left as they are. A rod therefore keeps its metals as
    separate nodes but collapses each linker to one vertex. This is the
    representation CrystalNets calls SingleNodes; MIL-53 gives bpq. An organic
    group that is itself periodic is left un-merged, so it cannot swallow a whole
    periodic chain.

    Arguments are the same as ``cgd_all_nodes``.
    '''
    edges, organic = _all_node_graph(
        atoms, components, breaking_pairs, regions, target_regions
    )
    edges = _merge_organic_nodes(edges, organic)
    return periodic_graph_cgd(dedup_periodic_edges(edges), name)


def _all_node_graph(atoms, components, breaking_pairs, regions, target_regions):
    '''
    Construct the all-node periodic graph.

    **returns:**
        edges: list of (u, v, sx, sy, sz) with 0-based node ids
        organic: dict node_id -> bool, True for carboxyl/linker (organic) nodes
                 and False for metal (inorganic) nodes
    '''
    tm = set(transition_metals())
    symbols = atoms.get_chemical_symbols()
    target_regions = set(int(r) for r in target_regions)

    kept_graph, kept_offsets = kept_bond_graph(atoms, breaking_pairs)
    self_translations = component_self_translations(
        components, kept_graph, kept_offsets
    )
    comp_images = component_atom_images(components, kept_graph, kept_offsets)

    atom_to_comp = {int(a): ci for ci, comp in enumerate(components) for a in comp}
    comp_to_region = {
        int(cid): int(rid) for rid, cids in regions.items() for cid in cids
    }
    target_comps = [
        c for c in range(len(components))
        if comp_to_region.get(c) in target_regions
    ]
    rod_comps = {c for c in target_comps if self_translations.get(c)}

    # atoms that carry a broken bond, grouped by their component
    broken_atoms = defaultdict(set)
    for pair in breaking_pairs:
        if len(pair) < 2:
            continue
        for atom in (int(pair[0]), int(pair[1])):
            broken_atoms[atom_to_comp[atom]].add(atom)

    # assign node ids: a discrete cluster is one node; a rod contributes one node
    # per metal atom and one per bridging (broken-bond) carbon
    node_of_key = {}
    organic = {}

    def node_id(key, is_organic):
        nid = node_of_key.setdefault(key, len(node_of_key))
        organic[nid] = is_organic
        return nid

    atom_node = {}
    for comp in target_comps:
        comp_atoms = set(int(a) for a in components[comp])
        if comp in rod_comps:
            for atom in comp_atoms:
                if symbols[atom] in tm:
                    atom_node[atom] = node_id(("metal", atom), False)
            for atom in broken_atoms[comp]:
                if symbols[atom] not in tm:
                    atom_node[atom] = node_id(("bridge", atom), True)
        else:
            cluster = node_id(("cluster", comp), False)
            for atom in comp_atoms:
                atom_node[atom] = cluster

    edges = []

    # within a rod, each non-node atom (the oxygens) contracts to edges between
    # the node atoms it bonds, using the real bond offsets so the periodic chain
    # bond survives
    for comp in rod_comps:
        comp_atoms = set(int(a) for a in components[comp])
        for oxygen in comp_atoms:
            if oxygen in atom_node:
                continue
            incidences = []
            for neigh in kept_graph.get(oxygen, []):
                neigh = int(neigh)
                if neigh not in atom_node:
                    continue
                for shift in kept_offsets.get((oxygen, neigh), []):
                    incidences.append((neigh, tuple(int(x) for x in shift)))
            for a in range(len(incidences)):
                for b in range(a + 1, len(incidences)):
                    (u_atom, su), (v_atom, sv) = incidences[a], incidences[b]
                    edges.append((
                        atom_node[u_atom], atom_node[v_atom],
                        su[0] - sv[0], su[1] - sv[1], su[2] - sv[2],
                    ))

    # linkers contract to an edge (ditopic) or a node (polytopic) between the
    # node atoms they attach to
    incidences = defaultdict(list)
    for pair in breaking_pairs:
        if len(pair) < 2:
            continue
        i, j = int(pair[0]), int(pair[1])
        shift = np.array([
            int(pair[2]) if len(pair) > 2 else 0,
            int(pair[3]) if len(pair) > 3 else 0,
            int(pair[4]) if len(pair) > 4 else 0,
        ], dtype=int)
        ci, cj = atom_to_comp[i], atom_to_comp[j]
        if i in atom_node and cj not in target_comps:
            node_atom, node_comp, link_atom, link_comp, sij = i, ci, j, cj, shift
        elif j in atom_node and ci not in target_comps:
            node_atom, node_comp, link_atom, link_comp, sij = j, cj, i, ci, -shift
        else:
            continue
        # a per-atom rod node is its own frame; a cluster node shares one frame
        origin = (np.zeros(3, dtype=int) if node_comp in rod_comps
                  else np.array(comp_images[node_comp].get(node_atom, (0, 0, 0))))
        link_image = np.array(comp_images[link_comp].get(link_atom, (0, 0, 0)))
        anchor = tuple((origin + sij - link_image).tolist())
        incidences[link_comp].append((atom_node[node_atom], anchor))

    next_linker_node = len(node_of_key)
    for link_comp, items in incidences.items():
        items = sorted(set(items))
        if len(items) < 2:
            continue
        if len(items) == 2:
            (u, su), (v, sv) = items
            edges.append((u, v, su[0] - sv[0], su[1] - sv[1], su[2] - sv[2]))
        else:
            lnode = next_linker_node
            next_linker_node += 1
            organic[lnode] = True
            for u, su in items:
                edges.append((u, lnode, su[0], su[1], su[2]))

    return edges, organic


def _merge_organic_nodes(edges, organic):
    '''
    Collapse each connected group of organic nodes into a single node.

    Metal (inorganic) nodes are untouched. A group of organic nodes joined by
    organic-organic edges becomes one node, with edge offsets adjusted for each
    member's image within the group. A group that is itself periodic (a bond
    closes back on a member in a different cell) is left alone, so a periodic
    organic chain is never swallowed into one point. This turns the all-node net
    into the single-node net.
    '''
    all_nodes = set()
    for u, v, *_ in edges:
        all_nodes.add(u)
        all_nodes.add(v)

    # adjacency among organic nodes via organic-organic edges, with offsets
    organic_adj = defaultdict(list)
    for u, v, sx, sy, sz in edges:
        if organic.get(u) and organic.get(v):
            organic_adj[u].append((v, (sx, sy, sz)))
            organic_adj[v].append((u, (-sx, -sy, -sz)))

    vmap = {}          # old node id -> new node id
    image_in_group = {}  # organic node -> its image within its merged group
    next_id = 0

    # inorganic nodes keep their own identity
    for node in sorted(all_nodes):
        if not organic.get(node):
            vmap[node] = next_id
            next_id += 1

    # organic nodes: group by connectivity, merge unless the group is periodic
    for seed in sorted(all_nodes):
        if not organic.get(seed) or seed in vmap:
            continue
        image = {seed: (0, 0, 0)}
        queue = [seed]
        periodic = False
        for node in queue:
            base = image[node]
            for neigh, off in organic_adj.get(node, []):
                new_off = (base[0] + off[0], base[1] + off[1], base[2] + off[2])
                if neigh in image:
                    if image[neigh] != new_off:
                        periodic = True
                else:
                    image[neigh] = new_off
                    queue.append(neigh)
        if periodic:
            for node in queue:
                vmap[node] = next_id
                next_id += 1
        else:
            group = next_id
            next_id += 1
            for node in queue:
                vmap[node] = group
                image_in_group[node] = image[node]

    remapped = []
    for u, v, sx, sy, sz in edges:
        du = image_in_group.get(u, (0, 0, 0))
        dv = image_in_group.get(v, (0, 0, 0))
        remapped.append((
            vmap[u], vmap[v],
            sx + du[0] - dv[0], sy + du[1] - dv[1], sz + du[2] - dv[2],
        ))
    return remapped


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
    ):
        '''
        Build a CGD PERIODIC_GRAPH using one of the supported methods.

        **parameters:**
            method: str
                One of:
                    - "sbus": each SBU is one node (rods collapse with a
                      self-edge).
                    - "all_node": rod SBUs keep their atoms as separate nodes,
                      recovering the true net (e.g. MIL-53 gives rna, not pcu).
                    - "single_node": the all-node net with each organic group
                      merged to one vertex (CrystalNets SingleNodes; MIL-53
                      gives bpq).

            name: str
                CGD graph ID.

        **returns:**
            cgd_text: str
                CGD PERIODIC_GRAPH content.
        '''
        guest_free = self.remove_guest()
        method = method.lower().strip()

        logger.info(
            "Building CGD with method='%s', name='%s'", method, name
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
                dedup_edges_for_systre=True,
            )

        elif method == "all_node":
            components, _, porphyrin, regions, breaking_pairs = mofdeconstructor.secondary_building_units(guest_free)

            if not components:
                raise RuntimeError("Deconstruction failed: no components (secondary_building_units).")
            if not regions:
                raise RuntimeError("Method 'all_node' requires regions, but regions is empty.")

            target_regions = regions_with_metal(guest_free, components, regions, porphyrin)
            if not target_regions:
                raise RuntimeError("No metal-containing regions found (excluding porphyrin metals) for method='all_node'.")

            return cgd_all_nodes(
                atoms=guest_free,
                components=components,
                breaking_pairs=breaking_pairs,
                regions=regions,
                target_regions=target_regions,
                name=name,
            )

        elif method == "single_node":
            components, _, porphyrin, regions, breaking_pairs = mofdeconstructor.secondary_building_units(guest_free)

            if not components:
                raise RuntimeError("Deconstruction failed: no components (secondary_building_units).")
            if not regions:
                raise RuntimeError("Method 'single_node' requires regions, but regions is empty.")

            target_regions = regions_with_metal(guest_free, components, regions, porphyrin)
            if not target_regions:
                raise RuntimeError("No metal-containing regions found (excluding porphyrin metals) for method='single_node'.")

            return cgd_single_nodes(
                atoms=guest_free,
                components=components,
                breaking_pairs=breaking_pairs,
                regions=regions,
                target_regions=target_regions,
                name=name,
            )

        else:
            raise ValueError("Unknown method. Choose from: 'sbus', 'all_node', 'single_node'.")

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
        choices=["sbus", "all_node", "single_node"],
        default="sbus",
        help="Topology extraction method."
    )
    parser.add_argument(
        "--name",
        default="net",
        help="CGD graph ID (default: net)."
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
        )
    except Exception as exc:
        logger.error("Failed to build CGD: %s", exc, exc_info=True)
        return 2

    TopologyExtractor.write_cgd(cgd_text, args.output)
    logger.info("CGD written to %s", args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())