#!/usr/bin/env python3
'''
Construction of CGD nets from crystal structures.

Systre works on nets rather than atoms, so a framework has to be reduced to
vertices and edges first. The node definition decides what the net describes,
and four are available: sbus places a vertex at each secondary building unit,
all_node keeps every branch point (splitting rod SBUs into their atoms), and
single_node coarsens all_node by merging each organic group to one vertex.
ligand_cluster represents complete organic ligands and metal clusters as the
two vertex classes of an incidence net. The same framework can give different
nets under each, which is expected rather than an error.
'''
from __future__ import annotations
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"

import argparse
import hashlib
import json
import logging
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
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


def _lattice_basis(vectors: Sequence[Sequence[int]]):
    '''
    Basis of the integer lattice generated by a set of translations.

    A rod running along c is discovered as (0, 0, 1) through one path and may
    surface as (0, 0, 2) through a longer one; a sheet yields two independent
    vectors mixed with their sums. What matters is the lattice they generate,
    not the individual vectors, so they are reduced to a row-echelon basis by
    integer (Hermite) elimination. The basis is canonical: the same lattice
    always gives the same rows whatever order the vectors arrive in.

    Dividing a single vector by its gcd would be wrong here -- a component that
    only closes on itself after two cells is periodic with (0, 0, 2), and
    calling that (0, 0, 1) invents a connection the structure does not have.

    **parameters:**
        vectors: list
            Integer translations (sx, sy, sz), duplicates allowed.

    **returns:**
        list
            Basis rows as tuples, ordered by leading non-zero column. Empty
            when the vectors generate nothing.
    '''
    rows = sorted(
        [int(v[0]), int(v[1]), int(v[2])] for v in vectors if any(int(x) for x in v)
    )
    basis = []

    for col in range(3):
        live = [k for k, row in enumerate(rows) if row[col] != 0]
        if not live:
            continue

        # Euclidean elimination: leave the gcd of this column in one row and
        # zero the column everywhere else.
        pivot = live[0]
        for other in live[1:]:
            while rows[other][col] != 0:
                if abs(rows[pivot][col]) > abs(rows[other][col]):
                    pivot, other = other, pivot
                factor = rows[other][col] // rows[pivot][col]
                rows[other] = [
                    rows[other][c] - factor * rows[pivot][c] for c in range(3)
                ]

        if rows[pivot][col] < 0:
            rows[pivot] = [-x for x in rows[pivot]]
        basis.append(rows[pivot])
        rows = [row for k, row in enumerate(rows) if k != pivot and any(row)]

    # reduce entries above each pivot so the basis is unique for the lattice
    for i in range(len(basis) - 1, -1, -1):
        pivot_col = next(c for c in range(3) if basis[i][c] != 0)
        for j in range(i):
            factor = basis[j][pivot_col] // basis[i][pivot_col]
            if factor:
                basis[j] = [basis[j][c] - factor * basis[i][c] for c in range(3)]

    return [tuple(row) for row in basis]


def _reduce_modulo_lattice(vec: Sequence[int], basis: Sequence[Tuple[int, int, int]]):
    '''
    Canonical representative of a translation modulo a lattice.

    A component that is periodic within itself is one vertex however far along
    its own chain a contact is made, so two contacts whose shifts differ by a
    translation of that component are the same edge of the net. Reducing both
    shifts to one representative is what makes the net independent of where the
    chain was unwrapped.

    **parameters:**
        vec: translation to reduce
        basis: lattice basis from ``_lattice_basis``

    **returns:**
        tuple
            (sx, sy, sz), unchanged when the basis is empty.
    '''
    reduced = [int(vec[0]), int(vec[1]), int(vec[2])]

    for row in basis:
        pivot_col = next(c for c in range(3) if row[c] != 0)
        factor = reduced[pivot_col] // row[pivot_col]
        if factor:
            reduced = [reduced[c] - factor * row[c] for c in range(3)]

    return (reduced[0], reduced[1], reduced[2])


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

    ``component_atom_images`` discards these back-edges; this recovers them. The
    translations found are reduced to a basis of the lattice they generate, so a
    rod yields one vector, a sheet two, and a discrete cluster none.

    **parameters:**
        components: list
            Connected components as lists of atom indices.

        kept_graph: python dictionary
            Atom-level graph after the broken bonds were removed.

        kept_offsets: python dictionary
            Periodic offsets of the kept bonds.

    **returns:**
        python dictionary
            Mapping of component id -> basis of the component's own translation
            lattice. Empty list for finite components.
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
                        if any(drift):
                            found.add(drift)

        self_translations[cid] = _lattice_basis(sorted(found))

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

    # A rod is one node however far along the chain a linker attaches, so the
    # shift to it is only defined modulo the rod's own translations. Unwrapping
    # a periodic component depends on where the traversal started, so without
    # this the same crystal written with its atoms in a different order gives a
    # different net -- MIL-53 came out pcu or msw depending on the ordering.
    pair_lattice = {}
    reduced_edges = []
    for u, v, sx, sy, sz in base_edges:
        if (u, v) not in pair_lattice:
            pair_lattice[(u, v)] = _lattice_basis(
                list(self_translations.get(u, [])) + list(self_translations.get(v, []))
            )
        reduced_edges.append(
            (u, v) + _reduce_modulo_lattice((sx, sy, sz), pair_lattice[(u, v)])
        )
    base_edges = reduced_edges

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


def coordination_contacts(
    atoms: Atoms,
    components: Sequence[Component],
    breaking_pairs: Sequence[BreakingPair],
    atom_graph: Dict[int, List[int]],
    is_cluster: Sequence[bool],
):
    '''
    Select the cut bonds that are genuine ligand-to-cluster coordination.

    ``ligands_and_metal_clusters`` reports a breaking pair for every rule that
    fires, and some of those pairs name two atoms that are not bonded to each
    other: the carboxylate rule, for instance, records the carboxylate carbon
    against the metal its oxygen is bound to. Roughly half the pairs returned
    for a carboxylate framework are of that kind, and they carry no real bond
    offset, so an edge built from one would sit at an arbitrary image. Only
    pairs that are bonds in the atom graph, and that join a ligand component to
    a cluster component, describe a contact.

    **parameters:**
        atoms: ASE atoms object

        components: connected components from the deconstructor

        breaking_pairs: broken bonds from the deconstructor

        atom_graph: atom-level connectivity of the intact structure

        is_cluster: per component, True when the component is a metal cluster

    **returns:**
        contacts: list of (metal_atom, donor_atom, cluster_id, ligand_id)

        atom_component: mapping atom index -> component id
    '''
    atom_component = {}
    for cid, comp in enumerate(components):
        for a in comp:
            atom_component[int(a)] = cid

    contacts = []
    seen = set()
    for pair in breaking_pairs:
        if len(pair) < 2:
            continue
        i, j = int(pair[0]), int(pair[1])
        if (i, j) in seen:
            continue
        if j not in [int(x) for x in atom_graph.get(i, [])]:
            continue

        u, v = atom_component.get(i), atom_component.get(j)
        if u is None or v is None or u == v:
            continue
        if is_cluster[u] == is_cluster[v]:
            continue

        seen.add((i, j))
        seen.add((j, i))
        if is_cluster[u]:
            contacts.append((i, j, u, v))
        else:
            contacts.append((j, i, v, u))

    return contacts, atom_component


def ligand_cluster_incidences(atoms: Atoms):
    '''
    Deconstruct a framework into ligand and cluster vertices and find every
    coordination incidence between them.

    The vertices are the components of ``ligands_and_metal_clusters``, which
    cuts the structure at the coordination bonds themselves: a ligand keeps
    every atom it owns, including its donor atoms, and a cluster keeps its
    metals with whatever inorganic bridges hold them together, so every atom
    belongs to exactly one vertex. An incidence is a contact between a ligand
    and one periodic image of a cluster, and it carries its denticity -- the
    number of donor bonds making that one contact -- because a carboxylate that
    has lost one oxygen to a defect still makes the same contact but with a
    lower denticity.

    A component that is periodic within itself, a rod or a sheet, is one vertex
    however far along its own chain a contact is made, so contact translations
    are reduced modulo that component's own translation lattice. Without this,
    chemically equivalent linkers of a rod framework come out with different
    coordination numbers and the answer changes when the same crystal is given
    in a supercell.

    **parameters:**
        atoms: ASE atoms object, guest-free

    **returns:**
        python dictionary with keys:
            components: the vertex components, as atom-index lists
            is_cluster: per component, True for a metal cluster
            lattices: per component, basis of its own translation lattice
            incidences: mapping (cluster, ligand, sx, sy, sz) -> denticity
    '''
    components, _, porphyrin, _, breaking_pairs = \
        mofdeconstructor.ligands_and_metal_clusters(atoms)
    if not components:
        raise ValueError("ligand_cluster: deconstruction returned no components.")

    atom_graph, _, bond_offsets = \
        mofdeconstructor.compute_ase_neighbour_with_offsets(atoms)

    tm = transition_metals()
    is_cluster = [
        component_has_transition_metal_not_in_porphyrin(atoms, comp, porphyrin, tm)
        for comp in components
    ]
    if not any(is_cluster):
        raise ValueError("ligand_cluster found no metal cluster.")
    if all(is_cluster):
        raise ValueError("ligand_cluster found no ligand coordinating a cluster.")

    contacts, _ = coordination_contacts(
        atoms, components, breaking_pairs, atom_graph, is_cluster
    )
    if not contacts:
        raise ValueError("ligand_cluster found no metal--ligand coordination contact.")

    # The kept graph is the framework with the coordination bonds cut, so each
    # component unwraps on its own bonds only.
    kept_graph, kept_offsets = remove_broken_bonds_from_breaking_pairs(
        atom_graph, bond_offsets, [(m, d) for m, d, _, _ in contacts]
    )
    images = component_atom_images(components, kept_graph, kept_offsets)
    lattices = component_self_translations(components, kept_graph, kept_offsets)

    pair_lattice = {}
    incidences = defaultdict(int)
    for metal, donor, cluster, ligand in contacts:
        if (cluster, ligand) not in pair_lattice:
            pair_lattice[(cluster, ligand)] = _lattice_basis(
                list(lattices.get(cluster, [])) + list(lattices.get(ligand, []))
            )

        shifts = bond_offsets.get((metal, donor))
        if not shifts:
            shifts = [
                (-s[0], -s[1], -s[2])
                for s in bond_offsets.get((donor, metal), [])
            ]
        if not shifts:
            shifts = [(0, 0, 0)]

        anchor = (
            np.array(images[cluster].get(metal, (0, 0, 0)), dtype=int)
            - np.array(images[ligand].get(donor, (0, 0, 0)), dtype=int)
        )
        for shift in shifts:
            translation = anchor + np.array(shift, dtype=int)
            key = (cluster, ligand) + _reduce_modulo_lattice(
                translation, pair_lattice[(cluster, ligand)]
            )
            incidences[key] += 1

    return {
        "components": components,
        "is_cluster": is_cluster,
        "lattices": lattices,
        "incidences": dict(incidences),
    }


def ligand_cluster_graph(atoms: Atoms, *, collapse_ditopic: bool = False):
    '''
    Construct the periodic incidence net between ligands and metal clusters.

    Vertices and edges come from ``ligand_cluster_incidences``: every complete
    ligand and every metal cluster is a vertex, and an edge records that a
    ligand coordinates one periodic image of a cluster. Several donor bonds
    within the same contact (a chelating or bridging carboxylate) are one edge,
    so chelation does not inflate the topological degree.

    A ligand held by a single contact -- coordinated solvent, a capping
    modulator, or a linker left dangling by a missing-linker defect -- is a
    vertex of degree one. No periodic net can carry one: Systre reports a
    collision and returns neither a name nor a relaxed net, so those ligands
    are left out of the graph. The defect is not lost, it shows in the reduced
    connectivity of the cluster that lost the linker, and
    ``ligand_cluster_fingerprint`` records the dangling ligand itself.

    **parameters:**
        atoms: ASE atoms object, guest-free

        collapse_ditopic: bool
            If True, a ligand bridging exactly two cluster contacts is spliced
            into a direct cluster--cluster edge instead of staying a vertex.
            The framework is unchanged, but the net is no longer subdivided, so
            it can carry an RCSR name (UiO-66 reads fcu rather than a 12,2-net
            that has no entry in the archive). Left False, the ligand stays a
            vertex and the net describes how each ligand meets the clusters.

    **returns:**
        edges: list of periodic graph edges

        ligand_nodes: set of node identifiers belonging to ligands

        node_atoms: mapping from node identifiers to represented atom indices
    '''
    deconstruction = ligand_cluster_incidences(atoms)
    components = deconstruction["components"]
    lattices = deconstruction["lattices"]
    incidences = set(deconstruction["incidences"])

    ligand_contacts = defaultdict(set)
    for cluster, ligand, sx, sy, sz in incidences:
        ligand_contacts[ligand].add((cluster, sx, sy, sz))
    terminal = {lig for lig, held in ligand_contacts.items() if len(held) < 2}
    if terminal:
        logger.info(
            "ligand_cluster: %d ligand(s) coordinate through a single contact "
            "(terminal solvent, modulator or dangling linker) and cannot be "
            "vertices of a net.",
            len(terminal),
        )
        incidences = {inc for inc in incidences if inc[1] not in terminal}
    if not incidences:
        raise ValueError(
            "ligand_cluster: every ligand is terminal, so there is no net."
        )

    participating_clusters = {cluster for cluster, _, *_ in incidences}
    participating_ligands = {ligand for _, ligand, *_ in incidences}
    cluster_map = {
        cluster: node for node, cluster in enumerate(sorted(participating_clusters))
    }
    ligand_start = len(cluster_map)
    ligand_map = {
        ligand: ligand_start + node
        for node, ligand in enumerate(sorted(participating_ligands))
    }

    edges = [
        (cluster_map[cluster], ligand_map[ligand], sx, sy, sz)
        for cluster, ligand, sx, sy, sz in incidences
    ]

    # give every periodic component its own connectivity back
    node_of_component = dict(cluster_map)
    node_of_component.update(ligand_map)
    for component, node in node_of_component.items():
        for sx, sy, sz in lattices.get(component, []):
            edges.append((node, node, sx, sy, sz))
    edges = dedup_periodic_edges(edges)

    ligand_nodes = {ligand_map[ligand] for ligand in participating_ligands}
    node_atoms = {
        node: set(int(a) for a in components[component])
        for component, node in node_of_component.items()
    }

    if collapse_ditopic:
        edges, spliced = splice_two_connected_nodes(edges, ligand_nodes)
        ligand_nodes -= spliced
        for node in spliced:
            node_atoms.pop(node, None)

    return edges, ligand_nodes, node_atoms


def splice_two_connected_nodes(edges, candidates):
    '''
    Replace each two-connected candidate vertex by a direct edge.

    A ditopic ligand is a connection rather than a branch point, and keeping it
    as a vertex subdivides the edge it makes. The subdivided net is a faithful
    description but has no RCSR name -- the archive holds pcu, not its
    subdivision -- so splicing restores the net that can be identified.

    **parameters:**
        edges: list of (u, v, sx, sy, sz)

        candidates: node ids allowed to be spliced out

    **returns:**
        edges: list of edges with the spliced vertices removed

        removed: set of node ids that were spliced out
    '''
    incident = defaultdict(list)
    self_looped = set()
    for u, v, sx, sy, sz in edges:
        if u == v:
            self_looped.add(u)
            continue
        incident[u].append((v, (sx, sy, sz)))
        incident[v].append((u, (-sx, -sy, -sz)))

    removed = {
        node for node in candidates
        if node not in self_looped and len(incident[node]) == 2
    }
    if not removed:
        return edges, removed

    kept = [
        (u, v, sx, sy, sz) for u, v, sx, sy, sz in edges
        if u not in removed and v not in removed
    ]
    for node in removed:
        (a, sa), (b, sb) = incident[node]
        kept.append((a, b, sa[0] - sb[0], sa[1] - sb[1], sa[2] - sb[2]))

    return dedup_periodic_edges(kept), removed


def cgd_ligand_cluster(atoms: Atoms, *, name: str = "net", collapse_ditopic: bool = False):
    '''Return a CGD representation of the ligand--metal-cluster incidence net.'''
    edges, _, _ = ligand_cluster_graph(atoms, collapse_ditopic=collapse_ditopic)
    return periodic_graph_cgd(edges, name)


def sbu_drawing_graph(atoms: Atoms):
    '''Construct the SBU-contracted graph together with its atom mapping.'''
    components, _, porphyrin, regions, breaking_pairs = \
        mofdeconstructor.secondary_building_units(atoms)
    target_regions = set(regions_with_metal(
        atoms, components, regions, porphyrin
    ))
    component_region = {
        int(component): int(region)
        for region, members in regions.items() for component in members
    }
    targets = [
        component for component in range(len(components))
        if component_region.get(component) in target_regions
    ]
    if not targets:
        raise ValueError("sbus found no metal-containing building unit.")

    kept_graph, kept_offsets = kept_bond_graph(atoms, breaking_pairs)
    base_edges = base_edges_with_shifts(
        atoms, components, breaking_pairs, kept_graph, kept_offsets
    )
    self_translations = component_self_translations(
        components, kept_graph, kept_offsets
    )
    pair_lattice = {}
    reduced_edges = []
    for u, v, sx, sy, sz in base_edges:
        if (u, v) not in pair_lattice:
            pair_lattice[(u, v)] = _lattice_basis(
                list(self_translations.get(u, []))
                + list(self_translations.get(v, []))
            )
        reduced_edges.append(
            (u, v) + _reduce_modulo_lattice(
                (sx, sy, sz), pair_lattice[(u, v)]
            )
        )
    base_edges = reduced_edges
    node_of_component = {
        component: node for node, component in enumerate(targets)
    }
    adjacency = defaultdict(list)
    for u, v, sx, sy, sz in base_edges:
        adjacency[u].append((v, (sx, sy, sz)))
        adjacency[v].append((u, (-sx, -sy, -sz)))

    edges = []
    organic = {node: False for node in node_of_component.values()}
    node_atoms = {
        node_of_component[component]: set(int(a) for a in components[component])
        for component in targets
    }
    next_node = len(node_of_component)
    for linker in range(len(components)):
        if linker in node_of_component:
            continue
        incidences = sorted(set(
            (node_of_component[neighbour], tuple(int(x) for x in shift))
            for neighbour, shift in adjacency.get(linker, [])
            if neighbour in node_of_component
        ))
        if len(incidences) < 2:
            continue
        if len(incidences) == 2:
            (u, su), (v, sv) = incidences
            edges.append((
                u, v, su[0] - sv[0], su[1] - sv[1], su[2] - sv[2]
            ))
            continue
        organic[next_node] = True
        node_atoms[next_node] = set(int(a) for a in components[linker])
        for u, shift in incidences:
            edges.append((u, next_node, *shift))
        next_node += 1

    for component, node in node_of_component.items():
        for shift in self_translations.get(component, []):
            edges.append((node, node, *shift))
    return dedup_periodic_edges(edges), organic, node_atoms


def _component_formula(atoms: Atoms, component: Component, divisor: int = 1):
    '''
    Hill-ordered formula of a component, used as its species label.

    A divisor greater than one writes the formula of a repeat unit rather than
    of the whole component, which is what a rod or a sheet needs: how much of
    an infinite cluster falls inside the cell depends on the cell, so the
    unreduced formula would change under a supercell.
    '''
    counts = defaultdict(int)
    for index in component:
        counts[atoms[int(index)].symbol] += 1

    ordered = []
    for symbol in ("C", "H"):
        if symbol in counts:
            ordered.append((symbol, counts.pop(symbol) // divisor))
    ordered.extend((symbol, n // divisor) for symbol, n in sorted(counts.items()))

    return "".join(
        symbol if n == 1 else f"{symbol}{n}" for symbol, n in ordered if n
    )


def _repeat_units(atoms: Atoms, component: Component, is_periodic: bool, n_contacts: int):
    '''
    How many repeat units of a component the cell holds.

    A finite cluster is one unit whatever the cell. A rod or a sheet is
    infinite, so the cell holds as many units as it holds: doubling the cell
    doubles its atoms and its contacts. The count is the largest divisor that
    divides the composition and the contacts alike, which recovers the repeat
    unit and makes everything derived from it independent of the cell.
    '''
    if not is_periodic:
        return 1

    counts = defaultdict(int)
    for index in component:
        counts[atoms[int(index)].symbol] += 1

    divisor = int(np.gcd.reduce(list(counts.values()) + [n_contacts]))
    return max(divisor, 1)


def _refine_colours(vertices, neighbours, initial, units, rounds: int = 3):
    '''
    Weisfeiler--Leman colour refinement on the ligand--cluster quotient graph.

    Counting how many ligands are 2-connected and how many clusters are
    12-connected does not say how those vertices sit relative to each other.
    Refinement does: each vertex repeatedly absorbs the colours of its
    neighbours, so after a few rounds a colour encodes the vertex's
    surroundings and the histogram of colours separates frameworks that share
    the same connectivity counts. Periodic shifts are deliberately ignored, and
    a neighbour count is taken per repeat unit, so the result does not depend
    on the cell the structure was written in.
    '''
    colours = dict(initial)

    for _ in range(rounds):
        refined = {}
        for vertex in vertices:
            histogram = Counter(colours[other] for other in neighbours[vertex])
            divisor = units.get(vertex, 1)
            refined[vertex] = (
                colours[vertex],
                tuple(sorted(
                    (colour, count // divisor if count % divisor == 0 else count)
                    for colour, count in histogram.items()
                )),
            )
        colours = refined

    return colours


def ligand_cluster_fingerprint(atoms: Atoms, *, algorithm: str = "sha256"):
    '''
    Describe how ligands meet metal clusters, without going through Systre.

    Systre answers a narrower question than this one. It names the net when the
    net has a name, but a ligand--cluster net keeps ditopic ligands as vertices
    and so is a subdivided net that RCSR does not list, and a framework with a
    dangling ligand is a net Systre refuses outright -- in both cases its
    topology hash is unavailable or, for an unnamed net, tied to the node
    numbering Systre happened to use. This fingerprint is computed from the
    deconstruction itself, so it exists for every framework, and it is
    invariant to atom order, to the choice of cell origin, and to being given a
    supercell of the same crystal.

    What it records is chemistry that a net alone cannot carry. Each ligand and
    cluster species is counted by formula, with a histogram of how many
    clusters each ligand bridges and of the denticity of those contacts. A
    missing-linker defect shows as a cluster species that has lost contacts; a
    linker left hanging by one end shows under ``terminal`` with its own
    formula, which is what separates it from a coordinated solvent; a
    carboxylate reduced from bridging to monodentate shows in the denticity
    histogram while the net itself is unchanged.

    **parameters:**
        atoms: ASE atoms object, guest-free

        algorithm: str
            Hash algorithm supported by hashlib.

    Every count is quoted per metal-cluster repeat unit, as an exact fraction:
    pristine HKUST-1 reads 4/3 of a C9H3O6 ligand per Cu2 paddlewheel, and one
    coordinated methanol in that cell reads 1/24.

    **returns:**
        python dictionary with keys:
            clusters: species -> {count, contacts histogram, capped histogram}
            ligands: species -> {count, contacts histogram, denticity histogram}
            terminal: species -> {count, denticity histogram} for ligands held
                by a single contact
            refinement: histogram of refined colours, as a sorted list of
                (colour digest, count)
            cluster_units: metal-cluster repeat units in the cell, the
                denominator the counts are quoted against
            fingerprint_hash: hex digest of everything above
    '''
    deconstruction = ligand_cluster_incidences(atoms)
    components = deconstruction["components"]
    is_cluster = deconstruction["is_cluster"]
    lattices = deconstruction["lattices"]
    incidences = deconstruction["incidences"]

    contacts = defaultdict(list)
    for (cluster, ligand, sx, sy, sz), denticity in incidences.items():
        contacts[cluster].append(((ligand, sx, sy, sz), denticity))
        contacts[ligand].append(((cluster, sx, sy, sz), denticity))

    vertices = sorted(contacts)

    # A ligand held by one contact caps a cluster rather than connecting it to
    # anything, so it is counted on its own and left out of the cluster's
    # framework connectivity. Coordinated solvent then does not change how
    # connected the cluster reads, while still being recorded.
    terminal_vertices = {
        vertex for vertex in vertices
        if not is_cluster[vertex] and len(contacts[vertex]) < 2
    }

    framework_vertices = [v for v in vertices if v not in terminal_vertices]
    neighbours = {
        vertex: [
            other for (other, *_), _ in contacts[vertex]
            if other not in terminal_vertices
        ]
        for vertex in framework_vertices
    }
    units = {
        vertex: _repeat_units(
            atoms, components[vertex], bool(lattices.get(vertex)), len(neighbours[vertex])
        )
        for vertex in framework_vertices
    }
    formula = {
        vertex: _component_formula(atoms, components[vertex], units.get(vertex, 1))
        for vertex in vertices
    }
    initial = {
        vertex: (
            "cluster" if is_cluster[vertex] else "ligand",
            formula[vertex],
            len(lattices.get(vertex, [])),
            len(neighbours[vertex]) // units[vertex],
            tuple(sorted(
                (denticity, count // units[vertex] if count % units[vertex] == 0 else count)
                for denticity, count in Counter(
                    denticity for (other, *_), denticity in contacts[vertex]
                    if other not in terminal_vertices
                ).items()
            )),
        )
        for vertex in framework_vertices
    }
    colours = _refine_colours(framework_vertices, neighbours, initial, units)

    clusters, ligands, terminal = {}, {}, {}
    for vertex in vertices:
        species = formula[vertex]
        divisor = units.get(vertex, 1)
        held_by = [
            (other, denticity) for (other, *_), denticity in contacts[vertex]
        ]
        framework = [
            (other, denticity) for other, denticity in held_by
            if other not in terminal_vertices
        ]

        if is_cluster[vertex]:
            entry = clusters.setdefault(
                species,
                {"count": 0, "contacts": defaultdict(int), "capped": defaultdict(int)},
            )
            entry["count"] += divisor
            entry["contacts"][len(framework) // divisor] += divisor
            entry["capped"][(len(held_by) - len(framework)) // divisor] += divisor
            continue

        if vertex in terminal_vertices:
            entry = terminal.setdefault(
                species, {"count": 0, "denticity": defaultdict(int)}
            )
        else:
            entry = ligands.setdefault(
                species,
                {"count": 0, "contacts": defaultdict(int), "denticity": defaultdict(int)},
            )
            entry["contacts"][len(held_by) // divisor] += divisor

        entry["count"] += divisor
        for _, denticity in held_by:
            entry["denticity"][denticity] += 1

    # weighted by repeat units, so one rod vertex holding two units of chain
    # counts as two, exactly as the two ligands bound to it do
    colour_counts = defaultdict(int)
    for vertex, colour in colours.items():
        colour_counts[colour] += units.get(vertex, 1)
    refinement = sorted(
        (_hash_text(repr(colour), algorithm=algorithm)[:16], count)
        for colour, count in colour_counts.items()
    )

    # Everything is quoted per metal-cluster repeat unit. A bigger box holds
    # proportionally more of everything, so the ratios do not move, while a
    # defect does move them. They are kept as exact fractions because one
    # capping solvent among twenty-four clusters must not round away to nothing,
    # which is what an integer count divided by a common factor would do.
    cluster_units = sum(
        units[vertex] for vertex in framework_vertices if is_cluster[vertex]
    )
    cluster_units = max(cluster_units, 1)

    def ratio(count):
        return str(Fraction(int(count), cluster_units))

    def scaled(histogram):
        return {int(k): ratio(v) for k, v in sorted(histogram.items())}

    def scale(group):
        return {
            species: {
                key: (ratio(value) if key == "count" else scaled(value))
                for key, value in entry.items()
            }
            for species, entry in sorted(group.items())
        }

    fingerprint = {
        "clusters": scale(clusters),
        "ligands": scale(ligands),
        "terminal": scale(terminal),
        "refinement": [[digest, ratio(count)] for digest, count in refinement],
    }
    # cluster_units says how big the cell was, not what is in it, so it is
    # reported but kept out of the hash: a supercell must hash to its crystal.
    fingerprint["fingerprint_hash"] = _hash_text(
        json.dumps(fingerprint, sort_keys=True, separators=(",", ":")),
        algorithm=algorithm,
    )
    fingerprint["cluster_units"] = cluster_units
    return fingerprint


def _hash_text(text: str, algorithm: str = "sha256"):
    '''Hex digest of a string.'''
    digest = hashlib.new(algorithm)
    digest.update(text.encode("utf-8"))
    return digest.hexdigest()


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
    edges, _, _ = _all_node_graph(
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
    edges, organic, _ = _all_node_graph(
        atoms, components, breaking_pairs, regions, target_regions
    )
    edges, _ = _merge_organic_nodes(edges, organic)
    return periodic_graph_cgd(dedup_periodic_edges(edges), name)


def _all_node_graph(atoms, components, breaking_pairs, regions, target_regions):
    '''
    Construct the all-node periodic graph.

    **returns:**
        edges: list of (u, v, sx, sy, sz) with 0-based node ids
        organic: dict node_id -> bool, True for carboxyl/linker (organic) nodes
                 and False for metal (inorganic) nodes
        node_atoms: dict node_id -> set of atom indices the node represents,
                 used to place the node at its real position when drawing
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
    node_atoms = defaultdict(set)
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
    for atom, nid in atom_node.items():
        node_atoms[nid].add(atom)

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
            node_atoms[lnode] = set(int(a) for a in components[link_comp])
            for u, su in items:
                edges.append((u, lnode, su[0], su[1], su[2]))

    return edges, organic, dict(node_atoms)


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
    return remapped, vmap


def net_geometry(atoms, *, method="all_node"):
    '''
    Trace the topological net over the real structure geometry.

    Each net node is placed at the real (Cartesian) centroid of the atoms it
    represents, so the returned net follows the actual framework rather than the
    idealised RCSR cell. Intended for drawing, not for topology identification.

    **parameters:**
        atoms: guest-free ASE atoms object

        method: "sbus", "all_node", "single_node" or "ligand_cluster"

    **returns:**
        positions: dict node_id -> (x, y, z) Cartesian centroid
        kinds: dict node_id -> "metal" or "organic"
        edges: list of (u, v, sx, sy, sz) with periodic shifts
        cell: 3x3 cell matrix (rows are lattice vectors)
    '''
    method = method.lower().strip()
    if method not in ("sbus", "all_node", "single_node", "ligand_cluster"):
        raise ValueError(
            "net_geometry supports 'sbus', 'all_node', 'single_node' or "
            "'ligand_cluster'."
        )

    if method == "ligand_cluster":
        edges, ligand_nodes, node_atoms = ligand_cluster_graph(atoms)
        organic = {node: node in ligand_nodes for node in node_atoms}
    elif method == "sbus":
        edges, organic, node_atoms = sbu_drawing_graph(atoms)
    else:
        components, _, porphyrin, regions, breaking_pairs = \
            mofdeconstructor.secondary_building_units(atoms)
        if not components:
            raise RuntimeError("Deconstruction produced no components.")
        target_regions = regions_with_metal(atoms, components, regions, porphyrin)
        if not target_regions:
            raise RuntimeError("No metal-containing regions found.")

        edges, organic, node_atoms = _all_node_graph(
            atoms, components, breaking_pairs, regions, target_regions
        )

    if method == "single_node":
        edges, vmap = _merge_organic_nodes(edges, organic)
        edges = dedup_periodic_edges(edges)
        merged_atoms = defaultdict(set)
        merged_organic = {}
        for nid, group in node_atoms.items():
            if nid not in vmap:
                continue
            merged_atoms[vmap[nid]].update(group)
            merged_organic[vmap[nid]] = (
                merged_organic.get(vmap[nid], False) or organic.get(nid, False)
            )
        node_atoms = dict(merged_atoms)
        organic = merged_organic

    cell = np.array(atoms.cell)
    frac = atoms.get_scaled_positions(wrap=False)

    positions = {}
    for nid, group in node_atoms.items():
        members = sorted(group)
        ref = frac[members[0]]
        unwrapped = []
        for atom in members:
            delta = frac[atom] - ref
            delta -= np.round(delta)  # minimum image relative to the first atom
            unwrapped.append(ref + delta)
        center = np.mean(unwrapped, axis=0)
        center -= np.floor(center)
        positions[nid] = np.asarray(center @ cell)

    kinds = {
        nid: ("organic" if organic.get(nid) else "metal") for nid in positions
    }
    edges = drawing_periodic_edges(positions, edges, cell)
    return positions, kinds, edges, cell


def drawing_periodic_edges(positions, edges, cell):
    '''
    Express periodic edges in the coordinate gauge used by drawn node
    positions.

    Periodic graph translations depend on the unit-cell representative chosen
    for every node. Node centroids are wrapped independently for display, so
    their representatives can differ from those used during graph
    construction. Integer gauge shifts are propagated along a spanning forest
    to reconcile the two choices without changing cycle translations. This
    preserves periodic self-edges and distinct connections between different
    images of the same pair of nodes.
    '''
    cell = np.asarray(cell, dtype=float)
    inv_cell = np.linalg.inv(cell)
    fractional = {
        node: np.asarray(position, dtype=float) @ inv_cell
        for node, position in positions.items()
    }
    adjacency = defaultdict(list)
    for u, v, sx, sy, sz in edges:
        if u == v:
            continue
        shift = np.array((sx, sy, sz), dtype=int)
        nearest = -np.rint(fractional[v] - fractional[u]).astype(int)
        adjacency[u].append((v, nearest - shift))
        adjacency[v].append((u, shift - nearest))

    gauge = {}
    for root in sorted(positions):
        if root in gauge:
            continue
        gauge[root] = np.zeros(3, dtype=int)
        queue = [root]
        for node in queue:
            for neighbour, delta in adjacency.get(node, []):
                if neighbour in gauge:
                    continue
                gauge[neighbour] = gauge[node] + delta
                queue.append(neighbour)

    adjusted = []
    for u, v, sx, sy, sz in edges:
        shift = (
            np.array((sx, sy, sz), dtype=int) + gauge[v] - gauge[u]
        )
        adjusted.append((u, v, int(shift[0]), int(shift[1]), int(shift[2])))
    return dedup_periodic_edges(adjusted)


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
                    - "ligand_cluster": complete organic ligands and metal
                      clusters form the two vertex classes of an incidence net.

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

        elif method == "ligand_cluster":
            return cgd_ligand_cluster(guest_free, name=name)

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
            raise ValueError(
                "Unknown method. Choose from: 'sbus', 'ligand_cluster', "
                "'all_node', 'single_node'."
            )

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
        choices=["sbus", "ligand_cluster", "all_node", "single_node"],
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
