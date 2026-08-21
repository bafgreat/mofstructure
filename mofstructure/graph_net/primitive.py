#!/usr/bin/env python3
'''
Reduction of a net to its primitive cell.

A crystal is very often described in a cell larger than it needs to be: a
1x1x2 supercell of a framework has twice as many vertex orbits as the net
requires, and every one of them is a translate of an orbit in the smaller
cell. Nothing in the vector representation says so, and two representations of
the same net in different cells are not equal as labelled quotient graphs, so
this has to be detected and undone before a canonical key can be trusted.

Delgado-Friedrichs & O'Keeffe make the detection easy by working in the
barycentric placement. A combinatorial translation of the net is an
automorphism, so it acts on the equilibrium placement as an affine map; being
a translation, it acts as an actual translation of that placement. So if the
net admits an extra translation carrying vertex orbit 0 to vertex orbit w,
the translation vector can only be

    t = p(w) - p(0),                                                     (1)

and there are at most n - 1 candidates to test. A candidate is a genuine
translation when it permutes the orbits and carries the edge set to itself.
Writing the induced map as u -> u' with

    p(u') = p(u) + t - c(u),   c(u) in Z^d,                              (2)

the covering vertex (u, m) goes to (u', m + c(u)), so an edge (u, v, s) is
carried to (u', v', s + c(v) - c(u)) and the test is whether every such image
is again an edge.

The valid translations together with Z^d generate a finer lattice L'. The net
is then rewritten on L': orbits collapse onto one representative each, and
each edge shift is re-expressed through

    s' = s + d(v) - d(u),                                                (3)

where d(u) is the translation taking the representative of u to u itself.
Because s, d(u) and d(v) all lie in L', the new shift is integral in a basis
of L'.

**references:**
    Delgado-Friedrichs, O. & O'Keeffe, M. (2003). Acta Cryst. A59, 351-360.
        doi:10.1107/S0108767303012017  (section 4, "our first task therefore
        is to find additional combinatorial translations, which must occur as
        actual translations of the equilibrium form")
    Delgado-Friedrichs, O. (2005). Discrete Comput. Geom. 33, 67-81.
        doi:10.1007/s00454-004-1147-x  (Corollary 9: an automorphism
        commuting with every translation acts as a translation)
'''
from __future__ import annotations

from fractions import Fraction
from collections.abc import Sequence

from mofstructure.graph_net.barycentric import barycentric_placement
from mofstructure.graph_net.lattice import (
    coords_in_rational_basis,
    rational_hermite_basis,
)
from mofstructure.graph_net.periodic_graph import PeriodicGraph

RatVec = tuple[Fraction, ...]


def _integer_part(value: Fraction) -> int:
    '''
    Floor of a rational, used to split a placement difference into lattice
    and residue parts.

    **parameters:**
        - value: fractions.Fraction

    **returns:**
        int
    '''
    return value.numerator // value.denominator


def additional_translations(
    graph: PeriodicGraph,
    placement: Sequence[RatVec] | None = None,
) -> list[tuple[RatVec, dict[int, int], dict[int, tuple[int, ...]]]]:
    '''
    Translations of the net that are not translations of the given cell.

    Applies equation (1) to generate candidates and equation (2) to test them.
    The identity is not included, so an empty result means the supplied cell is
    already primitive.

    **parameters:**
        - graph: PeriodicGraph
            Net whose quotient graph is connected.

        - placement: sequence or None
            Barycentric placement; computed when omitted.

    **returns:**
        list
            One (translation, orbit map, lattice correction) triple per
            non-trivial translation found.
    '''
    dim = graph.dim
    placement = placement if placement is not None else barycentric_placement(graph)
    edges = set(graph.edges)

    by_residue: dict[RatVec, int] = {}
    for v in range(graph.n_vertices):
        residue = tuple(x - _integer_part(x) for x in placement[v])
        by_residue.setdefault(residue, v)

    found = []
    for target in range(1, graph.n_vertices):
        shift = tuple(placement[target][k] - placement[0][k] for k in range(dim))

        mapping: dict[int, int] = {}
        correction: dict[int, tuple[int, ...]] = {}
        ok = True
        for u in range(graph.n_vertices):
            moved = tuple(placement[u][k] + shift[k] for k in range(dim))
            residue = tuple(x - _integer_part(x) for x in moved)
            image = by_residue.get(residue)
            if image is None:
                ok = False
                break
            mapping[u] = image
            correction[u] = tuple(
                _integer_part(moved[k] - placement[image][k] + Fraction(1, 2))
                for k in range(dim)
            )
        if not ok or len(set(mapping.values())) != graph.n_vertices:
            continue

        for (u, v, edge_shift) in graph.edges:
            new_shift = tuple(
                edge_shift[k] + correction[v][k] - correction[u][k] for k in range(dim)
            )
            candidate = (mapping[u], mapping[v], new_shift)
            if u > v:
                candidate = (mapping[v], mapping[u], tuple(-x for x in new_shift))
            normalised = PeriodicGraph.build(dim, graph.n_vertices, [candidate]).edges[0]
            if normalised not in edges:
                ok = False
                break
        if ok:
            found.append((shift, mapping, correction))
    return found


def primitive(graph: PeriodicGraph) -> PeriodicGraph:
    '''
    Rewrite a net on its primitive cell.

    Detects the additional translations, extends the lattice by them and
    collapses each orbit of vertices onto a single representative, re-expressing
    every edge through equation (3). A net already given primitively is
    returned unchanged.

    **parameters:**
        - graph: PeriodicGraph
            Net whose quotient graph is connected.

    **returns:**
        PeriodicGraph
            Equivalent net with the smallest number of vertex orbits and the
            finest translation lattice.
    '''
    dim = graph.dim
    if dim == 0 or graph.n_vertices == 1:
        return graph

    placement = barycentric_placement(graph)
    extra = additional_translations(graph, placement)
    if not extra:
        return graph

    generators: list[RatVec] = [
        tuple(Fraction(int(i == k)) for k in range(dim)) for i in range(dim)
    ]
    generators.extend(shift for (shift, _m, _c) in extra)
    basis = rational_hermite_basis(generators, dim)
    if len(basis) != dim:
        return graph

    parent = list(range(graph.n_vertices))

    def root(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for (_shift, mapping, _correction) in extra:
        for u, image in mapping.items():
            a, b = root(u), root(image)
            if a != b:
                parent[max(a, b)] = min(a, b)

    representative = {v: root(v) for v in range(graph.n_vertices)}
    order = sorted({representative[v] for v in range(graph.n_vertices)})
    index = {rep: i for i, rep in enumerate(order)}

    offset: dict[int, RatVec] = {}
    for v in range(graph.n_vertices):
        rep = representative[v]
        offset[v] = tuple(placement[v][k] - placement[rep][k] for k in range(dim))

    new_edges = []
    for (u, v, shift) in graph.edges:
        combined = tuple(
            Fraction(shift[k]) + offset[v][k] - offset[u][k] for k in range(dim)
        )
        coords = coords_in_rational_basis(combined, basis)
        if coords is None:
            return graph
        new_edges.append((index[representative[u]], index[representative[v]], coords))

    try:
        return PeriodicGraph.build(
            dim=dim,
            n_vertices=len(order),
            edges=new_edges,
            name=graph.name,
        )
    except ValueError:
        return graph


__all__ = ["additional_translations", "primitive"]
