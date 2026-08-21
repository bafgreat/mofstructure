#!/usr/bin/env python3
'''
Exact barycentric (equilibrium) placement of a periodic graph.

A placement assigns a point p(v) in R^d to every vertex. It is *barycentric*,
equivalently *in equilibrium*, when every vertex sits at the centre of gravity
of its neighbours,

    p(v) = ( 1 / |N(v)| ) * sum_{w in N(v)} p(w),                        (1)

where the sum runs over neighbours counted with multiplicity and each
neighbour is taken at its own image, so that a neighbour reached across an
edge (v, w, s) contributes p(w) + s. Rearranged, (1) is the statement that all
forces balance,

    sum_{(w, s) incident at v} ( p(w) + s - p(v) ) = 0.                  (2)

Collecting (2) over all vertices gives a linear system in the graph Laplacian
L = D - A of the quotient multigraph,

    ( L p )(v) = sum_{(w, s) incident at v} s.                           (3)

L has a one-dimensional kernel per coordinate for a connected graph, which is
exactly the freedom to translate the whole net, so pinning p(v_0) = 0 makes
the solution unique. Delgado-Friedrichs proves existence and uniqueness by
observing that a barycentric placement is the unique critical point of the
energy

    E(p) = sum_{edge orbits (v,w)} d( p(v), p(w) )^2,                    (4)

which is a non-degenerate quadratic form once one vertex is fixed, and whose
critical point is therefore a minimum. Two consequences matter downstream:
the placement minimises the sum of squared edge lengths, so it is an excellent
starting geometry that avoids gratuitous entanglement, and it is independent
of the metric, so it depends on the topology alone and not on the unit cell
the caller supplied.

Everything here is computed in exact rational arithmetic. This is not
fastidiousness: Delgado-Friedrichs & O'Keeffe warn that equilibrium positions
of distinct vertices can be arbitrarily close without being identical, so a
floating point comparison cannot decide whether two vertices collide
(Acta Cryst. A59, 351-360, section 4).

**stability.**
A net is *stable* when the placement is injective on the infinite vertex set,
so that a vertex is identified by its position. Vertices (u, t1) and (v, t2)
of the cover sit at p(u) + t1 and p(v) + t2, so they coincide precisely when

    p(u) - p(v)  is an integer vector.                                   (5)

A net is *locally stable* when, for each vertex, the neighbour images are
pairwise distinct. Local stability is the weaker condition the canonical form
algorithm actually needs, because it only ever sorts the neighbours of one
vertex at a time. Unstable nets are not merely awkward: a non-trivial
automorphism can act as the identity on the placement, so the combinatorial
symmetry group is strictly larger than any crystallographic group, and the
canonical form algorithm has no basis to order neighbours by. Systre refuses
these, and the ladder nets of Delgado-Friedrichs, O'Keeffe & Treacy are the
systematic family of examples.

**references:**
    Delgado-Friedrichs, O. (2005). Discrete Comput. Geom. 33, 67-81.
        doi:10.1007/s00454-004-1147-x  (Theorem 4: existence and uniqueness)
    Delgado-Friedrichs, O. & O'Keeffe, M. (2003). Acta Cryst. A59, 351-360.
        doi:10.1107/S0108767303012017  (section 3, equations 1 and 2)
    Delgado-Friedrichs, O. (2004). Graph Drawing (GD 2003), LNCS 2912, 178-189.
        (section 3: stability and local stability)
    Delgado-Friedrichs, O., O'Keeffe, M. & Treacy, M. M. J. (2020).
        Acta Cryst. A76, 735-738. doi:10.1107/S2053273320012905  (ladder nets)
'''
from __future__ import annotations

from fractions import Fraction
from collections.abc import Sequence

from mofstructure.graph_net.periodic_graph import PeriodicGraph

RatVec = tuple[Fraction, ...]


class UnstableNetError(ValueError):
    '''
    Raised when a net is not locally stable and cannot be canonicalised.

    **parameters:**
        - message: str
            Human readable description.

        - collisions: list
            List of colliding items, either vertex pairs or (vertex, position)
            groups depending on which check failed.
    '''

    def __init__(self, message: str, collisions: list | None = None):
        super().__init__(message)
        self.collisions = collisions or []


def _solve_exact(matrix: list[list[Fraction]], rhs: list[list[Fraction]]) -> list[list[Fraction]]:
    '''
    Solve `matrix @ x = rhs` exactly by Gauss-Jordan elimination over Q.

    Several right hand sides are solved at once, one per Cartesian coordinate.

    **parameters:**
        - matrix: list
            Square coefficient matrix of Fraction, modified in place.

        - rhs: list
            Right hand side, one row per equation and one column per
            coordinate.

    **returns:**
        list
            Solution rows, one per unknown.

    **raises:**
        ValueError
            If the matrix is singular, which for a reduced Laplacian means the
            quotient graph was not connected.
    '''
    n = len(matrix)
    width = len(rhs[0]) if rhs else 0
    for col in range(n):
        pivot = None
        for row in range(col, n):
            if matrix[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            raise ValueError("singular system: the quotient graph is not connected")
        matrix[col], matrix[pivot] = matrix[pivot], matrix[col]
        rhs[col], rhs[pivot] = rhs[pivot], rhs[col]
        inv = Fraction(1, 1) / matrix[col][col]
        matrix[col] = [x * inv for x in matrix[col]]
        rhs[col] = [x * inv for x in rhs[col]]
        for row in range(n):
            if row == col:
                continue
            factor = matrix[row][col]
            if factor == 0:
                continue
            matrix[row] = [
                matrix[row][k] - factor * matrix[col][k] for k in range(n)
            ]
            rhs[row] = [rhs[row][k] - factor * rhs[col][k] for k in range(width)]
    return rhs


def barycentric_placement(graph: PeriodicGraph) -> list[RatVec]:
    '''
    Exact equilibrium placement with the first vertex pinned at the origin.

    Builds and solves the Laplacian system (3). Self-loops drop out of the
    left hand side, because a loop (v, v, s) contributes both +s and -s to the
    incidences of v and adds 2 to both the degree and the diagonal of A, so it
    can never influence the equilibrium position of its own vertex.

    **parameters:**
        - graph: PeriodicGraph
            Quotient graph with a connected quotient.

    **returns:**
        list
            One tuple of `graph.dim` Fractions per vertex orbit.

    **raises:**
        ValueError
            If the quotient graph is not connected.
    '''
    if not graph.is_quotient_connected():
        raise ValueError("barycentric placement requires a connected quotient graph")
    n = graph.n_vertices
    d = graph.dim
    zero = tuple(Fraction(0) for _ in range(d))
    if n == 1 or d == 0:
        return [zero for _ in range(n)]

    inc = graph.incidences()
    size = n - 1
    matrix = [[Fraction(0) for _ in range(size)] for _ in range(size)]
    rhs = [[Fraction(0) for _ in range(d)] for _ in range(size)]

    for v in range(1, n):
        row = v - 1
        for (w, shift) in inc[v]:
            matrix[row][row] += 1
            if w != 0:
                matrix[row][w - 1] -= 1
            for k in range(d):
                rhs[row][k] += shift[k]

    solution = _solve_exact(matrix, rhs)
    placement = [zero]
    for row in solution:
        placement.append(tuple(row))
    return placement


def edge_vectors(graph: PeriodicGraph, placement: Sequence[RatVec]) -> list[tuple[int, int, RatVec]]:
    '''
    Directed edge vectors of the barycentric drawing.

    For every edge orbit (v, w, s) both directions are emitted, the forward
    vector being

        e = p(w) + s - p(v)

    and the reverse being its negation. Algorithm 5 draws its candidate
    coordinate bases from exactly this list.

    **parameters:**
        - graph: PeriodicGraph
            Quotient graph.

        - placement: sequence
            Barycentric placement, as returned by `barycentric_placement`.

    **returns:**
        list
            List of (tail, head, vector) triples, both directions per edge.
    '''
    out = []
    for (u, v, shift) in graph.edges:
        vec = tuple(placement[v][k] + shift[k] - placement[u][k] for k in range(graph.dim))
        out.append((u, v, vec))
        out.append((v, u, tuple(-x for x in vec)))
    return out


def find_collisions(graph: PeriodicGraph, placement: Sequence[RatVec]) -> list[tuple[int, int]]:
    '''
    Vertex orbit pairs that occupy the same point of the barycentric drawing.

    Applies criterion (5): orbits u and v collide when every component of
    p(u) - p(v) is an integer.

    **parameters:**
        - graph: PeriodicGraph
            Quotient graph.

        - placement: sequence
            Barycentric placement.

    **returns:**
        list
            Sorted list of colliding (u, v) pairs with u < v.
    '''
    out = []
    n = graph.n_vertices
    for u in range(n):
        for v in range(u + 1, n):
            diff = [placement[u][k] - placement[v][k] for k in range(graph.dim)]
            if all(x.denominator == 1 for x in diff):
                out.append((u, v))
    return out


def find_local_collisions(
    graph: PeriodicGraph,
    placement: Sequence[RatVec],
) -> list[tuple[int, RatVec]]:
    '''
    Vertices whose neighbour images are not pairwise distinct.

    Local stability is what Algorithm 5 relies on when it sorts the
    neighbours of a vertex into lexicographic order: if two neighbours share a
    position the order is arbitrary and the resulting key is not canonical.

    **parameters:**
        - graph: PeriodicGraph
            Quotient graph.

        - placement: sequence
            Barycentric placement.

    **returns:**
        list
            List of (vertex, duplicated neighbour position) pairs.
    '''
    out = []
    inc = graph.incidences()
    for v in range(graph.n_vertices):
        seen = {}
        for (w, shift) in inc[v]:
            pos = tuple(placement[w][k] + shift[k] for k in range(graph.dim))
            if pos in seen:
                out.append((v, pos))
            seen[pos] = w
    return out


def is_stable(graph: PeriodicGraph, placement: Sequence[RatVec] | None = None) -> bool:
    '''
    Whether the barycentric placement is injective on the infinite net.

    **parameters:**
        - graph: PeriodicGraph
            Quotient graph.

        - placement: sequence or None
            Precomputed placement; recomputed when omitted.

    **returns:**
        bool
    '''
    placement = placement if placement is not None else barycentric_placement(graph)
    return not find_collisions(graph, placement)


def is_locally_stable(graph: PeriodicGraph, placement: Sequence[RatVec] | None = None) -> bool:
    '''
    Whether every vertex has pairwise distinct neighbour positions.

    **parameters:**
        - graph: PeriodicGraph
            Quotient graph.

        - placement: sequence or None
            Precomputed placement; recomputed when omitted.

    **returns:**
        bool
    '''
    placement = placement if placement is not None else barycentric_placement(graph)
    return not find_local_collisions(graph, placement)


__all__ = [
    "UnstableNetError",
    "barycentric_placement",
    "edge_vectors",
    "find_collisions",
    "find_local_collisions",
    "is_stable",
    "is_locally_stable",
]
