#!/usr/bin/env python3
'''
Canonical form of a periodic graph: the Systre key.

A net has infinitely many vector representations, because one may relabel the
vertex orbits, translate each representative to another of its images, and
apply any unimodular change of lattice basis. Deciding whether two nets are
the same therefore needs a *canonical form*: a representation that depends
only on the isomorphism type. Once two nets are reduced to canonical form,
testing them for equality is string comparison, and a database of nets becomes
a dictionary keyed on that string. This is exactly how the RCSR archive that
ships with Systre is organised.

The construction implemented here is Algorithm 5 of Delgado-Friedrichs,
*Barycentric Drawings of Periodic Graphs* (GD 2003, LNCS 2912, page 187),
which is the algorithm behind Systre. The idea is that a representation is
pinned down by two choices, a start vertex and a basis, and that the
barycentric placement supplies a canonical way to make every other choice:

  1. Compute the exact barycentric placement p. It is unique up to affine
     transformation and depends on the topology alone.

  2. Form the set A of directed edge vectors e = p(w) + s - p(v) of that
     drawing, both directions per edge orbit.

  3. For each ordered d-tuple of linearly independent vectors from A, let B be
     the matrix with those vectors as columns and traverse the graph breadth
     first from the tail of the first vector, visiting the neighbours of each
     vertex in lexicographic order of B^-1 e. Because the net is locally
     stable, that order is unambiguous, so the traversal is fully determined
     by the tuple.

  4. Each traversal numbers the vertices in visiting order and writes out one
     record (I(v), I(w), t) per edge. Keep the lexicographically smallest such
     list over all tuples.

  5. The accumulated shifts t live in the basis B, which is generally not a
     basis of the translation lattice, so finally re-express them in the
     canonical Hermite basis of the lattice they generate. That returns
     integers, and the result is the key.

Step 3 is quadratic-to-cubic in the number of edges if taken literally. The
practical reduction given by Delgado-Friedrichs & O'Keeffe is used here: the
tuple is required to start from a common vertex, all d of them where the
coordination allows it, and otherwise the longest prefix that does. Their own
example of a net that forces the fallback is nbo, in which every vertex is
square planar so no three edge vectors at a vertex are independent
(Acta Cryst. A59, 351-360, section 5).

**two deviations from the printed pseudocode.**
The pseudocode has been checked against all 2930 keys in the shipped RCSR
archive, and two conventions in it do not reproduce that archive:

  * the self-loop sign rule is printed as ``if v = w and t < 0: t = -t``,
    while the archive stores loops with a *negative* leading component. See
    `periodic_graph._normalise_edge` for the evidence.

  * the pseudocode selects the smallest traversal before re-expressing the
    shifts in an integral basis. Doing so leaves the key on whatever basis the
    winning tuple happened to induce, which is not canonical: two
    representations of dia that agree on the traversal can still differ by a
    unimodular map. The minimum is therefore taken over the *finished* key,
    after the integral basis has been imposed, which is the only choice that
    makes the result independent of the input representation.

**references:**
    Delgado-Friedrichs, O. (2004). Graph Drawing (GD 2003), LNCS 2912, 178-189.
        Springer.  (Algorithm 5, page 187)
    Delgado-Friedrichs, O. & O'Keeffe, M. (2003). Acta Cryst. A59, 351-360.
        doi:10.1107/S0108767303012017  (sections 4 and 5)
    Delgado-Friedrichs, O. (2005). Discrete Comput. Geom. 33, 67-81.
        doi:10.1007/s00454-004-1147-x
'''
from __future__ import annotations

from collections import deque
from fractions import Fraction
from itertools import permutations
from math import gcd
from collections.abc import Sequence

from mofstructure.graph_net.barycentric import (
    UnstableNetError,
    barycentric_placement,
    edge_vectors,
    find_local_collisions,
)
from mofstructure.graph_net.lattice import (
    coords_in_basis,
    hermite_basis,
    hermite_basis_with_transform,
)
from mofstructure.graph_net.periodic_graph import PeriodicGraph
from mofstructure.graph_net.primitive import primitive

IntVec = tuple[int, ...]
Record = tuple[int, int, IntVec]


def _rank(vectors: Sequence[Sequence[int]], dim: int) -> int:
    '''
    Rank of a set of integer vectors, by exact elimination over Q.

    **parameters:**
        - vectors: sequence
            Sequence of integer vectors.

        - dim: int
            Ambient dimension.

    **returns:**
        int
    '''
    rows = [[Fraction(x) for x in v] for v in vectors]
    rank = 0
    for col in range(dim):
        pivot = None
        for i in range(rank, len(rows)):
            if rows[i][col] != 0:
                pivot = i
                break
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inv = Fraction(1) / rows[rank][col]
        rows[rank] = [x * inv for x in rows[rank]]
        for i in range(len(rows)):
            if i != rank and rows[i][col] != 0:
                factor = rows[i][col]
                rows[i] = [rows[i][k] - factor * rows[rank][k] for k in range(dim)]
        rank += 1
        if rank == len(rows):
            break
    return rank


def _det(matrix: Sequence[Sequence[int]], dim: int) -> int:
    '''
    Determinant of a small integer matrix by cofactor expansion.

    **parameters:**
        - matrix: sequence
            Rows of the matrix.

        - dim: int
            Side length.

    **returns:**
        int
    '''
    if dim == 1:
        return matrix[0][0]
    if dim == 2:
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    total = 0
    for c in range(dim):
        if matrix[0][c] == 0:
            continue
        minor = [
            [matrix[r][k] for k in range(dim) if k != c]
            for r in range(1, dim)
        ]
        total += ((-1) ** c) * matrix[0][c] * _det(minor, dim - 1)
    return total


def _scaled_inverse(columns: Sequence[IntVec], dim: int) -> tuple[list[list[int]], int] | None:
    '''
    Adjugate of the matrix whose columns are `columns`, with its determinant.

    Working with the adjugate keeps the traversal in exact integer arithmetic.
    Since

        B^-1 = adj(B) / det(B),

    applying `sign(det) * adj(B)` to an edge vector gives `|det| * B^-1 e`,
    which is a positive multiple of the true coordinate vector and therefore
    induces exactly the same lexicographic order on the neighbours of a
    vertex. The uniform factor |det| cancels later, because the shifts are
    re-expressed in a basis drawn from the shifts themselves.

    **parameters:**
        - columns: sequence
            `dim` integer vectors taken as the columns of the matrix.

        - dim: int
            Side length.

    **returns:**
        tuple or None
            (rows of sign(det) * adj(B), det), or None when singular.
    '''
    mat = [[columns[c][r] for c in range(dim)] for r in range(dim)]
    det = _det(mat, dim)
    if det == 0:
        return None
    sign = 1 if det > 0 else -1
    adj = [[0] * dim for _ in range(dim)]
    for r in range(dim):
        for c in range(dim):
            minor = [
                [mat[rr][cc] for cc in range(dim) if cc != r]
                for rr in range(dim) if rr != c
            ]
            cof = ((-1) ** (r + c)) * (_det(minor, dim - 1) if dim > 1 else 1)
            adj[r][c] = sign * cof
    return adj, det


def _apply(matrix: Sequence[Sequence[int]], vec: Sequence[int], dim: int) -> IntVec:
    '''
    Multiply an integer matrix by an integer column vector.

    **parameters:**
        - matrix: sequence
            Rows of the matrix.

        - vec: sequence
            Column vector.

        - dim: int
            Dimension.

    **returns:**
        tuple
    '''
    if dim == 3:
        row0, row1, row2 = matrix
        v_0, v_1, v_2 = vec
        return (
            row0[0] * v_0 + row0[1] * v_1 + row0[2] * v_2,
            row1[0] * v_0 + row1[1] * v_1 + row1[2] * v_2,
            row2[0] * v_0 + row2[1] * v_1 + row2[2] * v_2,
        )
    if dim == 2:
        row0, row1 = matrix
        v_0, v_1 = vec
        return (row0[0] * v_0 + row0[1] * v_1, row1[0] * v_0 + row1[1] * v_1)
    return tuple(
        sum(matrix[r][c] * vec[c] for c in range(dim)) for r in range(dim)
    )


def _candidate_bases(
    directed: Sequence[tuple[int, int, IntVec]],
    dim: int,
) -> list[tuple[int, tuple[IntVec, ...]]]:
    '''
    Enumerate the ordered tuples of edge vectors used as coordinate bases.

    Following the reduction of Delgado-Friedrichs & O'Keeffe (2003, section 5)
    the tuple is anchored at a vertex: the longest possible prefix is drawn
    from the directed edges leaving one common vertex, and only the remainder,
    when the coordination figure is degenerate, is drawn from the whole edge
    set. Their own example of a net forcing the fallback is nbo, whose
    vertices are square planar so that no three edge vectors at a vertex are
    independent.

    **parameters:**
        - directed: sequence
            Directed edge vectors as (tail, head, vector).

        - dim: int
            Dimension of the net.

    **returns:**
        list
            List of (start_vertex, tuple_of_vectors) candidates.
    '''
    at_vertex: dict[int, list[IntVec]] = {}
    for (tail, _head, vec) in directed:
        at_vertex.setdefault(tail, []).append(vec)

    best_prefix = max(
        (min(_rank(vecs, dim), dim) for vecs in at_vertex.values()),
        default=0,
    )
    if best_prefix == 0:
        return []

    all_vectors = [vec for (_t, _h, vec) in directed]
    candidates: list[tuple[int, tuple[IntVec, ...]]] = []
    for start, vecs in sorted(at_vertex.items()):
        for prefix in permutations(vecs, best_prefix):
            if _rank(prefix, dim) < best_prefix:
                continue
            if best_prefix == dim:
                candidates.append((start, tuple(prefix)))
                continue
            for suffix in permutations(all_vectors, dim - best_prefix):
                whole = tuple(prefix) + tuple(suffix)
                if _rank(whole, dim) == dim:
                    candidates.append((start, whole))
    return candidates


def _traverse(
    start: int,
    transform: Sequence[Sequence[int]],
    incidences: dict[int, list[tuple[int, IntVec, IntVec]]],
    n_vertices: int,
    dim: int,
) -> tuple[list[tuple[int, int, IntVec, IntVec]], dict[int, int]] | None:
    '''
    Breadth-first traversal producing one candidate representation.

    Vertices are numbered 1, 2, ... in visiting order. At each vertex the
    incident edges are taken in lexicographic order of their edge vector
    expressed in the chosen basis, which is well defined precisely because the
    net is locally stable. For each incident edge the accumulated translation

        t = q(v) + e - q(w)

    is recorded together with the pair of visiting numbers, once per edge: the
    guard keeps the record when the head has the larger number, and for a loop
    keeps the direction with positive t.

    **parameters:**
        - start: int
            Vertex to start from.

        - transform: sequence
            Rows of sign(det) * adj(B) for the chosen basis B.

        - incidences: python dictionary
            Vertex -> list of (neighbour, edge vector).

        - n_vertices: int
            Expected number of vertex orbits.

        - dim: int
            Dimension.

    **returns:**
        tuple or None
            (records, numbering) where records are (number, number,
            translation) triples and numbering maps each vertex orbit to its
            visiting number, or None when the traversal failed to reach every
            vertex.
    '''
    zero = (0,) * dim
    number = {start: 1}
    potential = {start: zero}
    lattice_potential = {start: zero}
    next_number = 2
    queue = deque([start])
    records: list[tuple[int, int, IntVec, IntVec]] = []

    while queue:
        v = queue.popleft()
        transformed = [
            (w, _apply(transform, vec, dim), raw) for (w, vec, raw) in incidences[v]
        ]
        transformed.sort(key=lambda item: item[1])
        for (w, vec, raw) in transformed:
            if w not in number:
                number[w] = next_number
                potential[w] = tuple(potential[v][k] + vec[k] for k in range(dim))
                lattice_potential[w] = tuple(
                    lattice_potential[v][k] + raw[k] for k in range(dim)
                )
                next_number += 1
                queue.append(w)
            shift = tuple(
                potential[v][k] + vec[k] - potential[w][k] for k in range(dim)
            )
            lattice_shift = tuple(
                lattice_potential[v][k] + raw[k] - lattice_potential[w][k]
                for k in range(dim)
            )
            if number[v] < number[w] or (number[v] == number[w] and shift > zero):
                records.append((number[v], number[w], shift, lattice_shift))

    if len(number) != n_vertices:
        return None
    return records, number


def _coords(basis: Sequence[IntVec], vec: Sequence[int], dim: int) -> IntVec | None:
    '''
    Express `vec` as an integer combination of `basis`.

    Solves x . basis = vec over Q and rejects non-integral solutions, which is
    how a basis generating only a proper sublattice is detected.

    **parameters:**
        - basis: sequence
            `dim` integer vectors forming a basis, as rows.

        - vec: sequence
            Vector to re-express.

        - dim: int
            Dimension.

    **returns:**
        tuple or None
    '''
    mat = [[Fraction(basis[r][c]) for r in range(dim)] for c in range(dim)]
    rhs = [Fraction(x) for x in vec]
    for col in range(dim):
        pivot = None
        for row in range(col, dim):
            if mat[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            return None
        mat[col], mat[pivot] = mat[pivot], mat[col]
        rhs[col], rhs[pivot] = rhs[pivot], rhs[col]
        scale = Fraction(1) / mat[col][col]
        mat[col] = [x * scale for x in mat[col]]
        rhs[col] = rhs[col] * scale
        for row in range(dim):
            if row == col or mat[row][col] == 0:
                continue
            factor = mat[row][col]
            mat[row] = [mat[row][k] - factor * mat[col][k] for k in range(dim)]
            rhs[row] = rhs[row] - factor * rhs[col]
    if any(x.denominator != 1 for x in rhs):
        return None
    return tuple(int(x) for x in rhs)


def _finalise_all(
    records: Sequence[tuple[int, int, IntVec, IntVec]],
    dim: int,
    best: Sequence[Record] | None = None,
) -> tuple[tuple[Record, ...], int] | None:
    '''
    Re-express a traversal over the integers.

    The translations produced by `_traverse` live in the coordinates induced
    by the traversal's own frame, which is built from edge vectors rather than
    lattice vectors, so they have to be re-expressed before they can be
    written as a key. Which integral basis is used is not free: the infinitely
    many unimodular matrices would otherwise give infinitely many equally
    valid strings for one net.

    The basis taken here is the Hermite normal form of the lattice in which
    the traversal's translations live. Being a function of the lattice alone it
    is canonical, and being a genuine basis of that lattice it always exists
    and always expresses every translation over the integers.

    An attempt was made to derive this lattice cheaply from the d columns of
    the traversal map instead of from the m record translations, on the
    argument that the records span the image of Z^d. Reducing d vectors rather
    than one per edge would turn an O(m d^2) step into an O(d^3) one, which is
    the dominant cost for a large net. It is recorded here only as a warning:
    the two lattices are *not* equal in general, the substitution changes keys,
    and it was reverted. Any future attempt needs to establish the relationship
    between the two lattices first, not assume it.

    An earlier version instead adopted the first `dim` linearly independent
    translations in traversal order, which is the rule Systre itself follows
    and which reproduces stored archive keys exactly. It was abandoned because
    it is not total: independence does not imply generation, so those `dim`
    vectors can span a proper sublattice, and for eighteen RCSR nets, among
    them cdj, css and elv, *every* candidate traversal ran into that and no
    key could be produced at all. A lattice of rank d need not have a basis
    among any d of its generators, so the failure is intrinsic to the rule
    rather than a matter of searching harder.

    **parameters:**
        - records: sequence
            Traversal records carrying both the transformed and the lattice
            translation, in traversal order.

        - dim: int
            Dimension.

    Only the basis itself is needed to compare one candidate against another,
    so the reduction is run without tracking the unimodular transform that
    produced it. That matters: carrying the transform widens every row of the
    reduction from `dim` to `dim + m` entries, and with one generator per edge
    it is the single most expensive step in the whole canonical form. The
    transform is needed only to express the basis back in the original lattice
    coordinates, which only the winning traversals require, so `canonical_form`
    recovers it from `finalisation_basis` once the winners are known.

        best: sequence or None
            The finished key of the incumbent, when there is one. Supplying it
            lets the routine stop as soon as this candidate has lost, which
            does not change the winner: the comparison is decided by the first
            record that differs and the records are produced in sorted order.

    **returns:**
        tuple or None
            (sorted integer records, ordering) where ordering is -1 when this
            candidate beat the incumbent and 0 when it tied. None when the
            translations fail to span a lattice of full rank, or when the
            candidate lost.
    '''
    shifts = [rec[2] for rec in records]
    basis = hermite_basis(shifts, dim)
    if len(basis) != dim:
        return None

    # The finished key is sorted on (a, b, coords), so its blocks of equal
    # (a, b) fall in an order the traversal already fixes, and coordinates
    # settle only the order inside a block. Walking the blocks in that
    # order and re-expressing one block at a time visits the records in
    # exactly the order the sorted key holds them, which is what allows a
    # candidate to be dropped at the first record that loses to the incumbent
    # without changing which candidate wins.
    order = sorted(range(len(records)), key=lambda i: records[i][:2])
    out: list[Record] = []
    verdict = 0
    head = 0
    while head < len(order):
        a, b = records[order[head]][:2]
        tail = head
        while tail < len(order) and records[order[tail]][:2] == (a, b):
            tail += 1
        block: list[Record] = []
        for position in order[head:tail]:
            shift = records[position][2]
            coords = coords_in_basis(shift, basis)
            if coords is None:
                return None
            if a == b:
                lead = next((x for x in coords if x != 0), 0)
                if lead > 0:
                    coords = tuple(-x for x in coords)
            block.append((a, b, coords))
        block.sort()
        for record in block:
            # Once the candidate is ahead the comparison is settled, but the
            # remaining records are still needed: they become the incumbent.
            if best is not None and verdict == 0:
                incumbent = best[len(out)]
                if record < incumbent:
                    verdict = -1
                elif record > incumbent:
                    return None
            out.append(record)
        head = tail
    return tuple(out), verdict


def finalisation_basis(
    records: Sequence[tuple[int, int, IntVec, IntVec]],
    dim: int,
) -> list[IntVec] | None:
    '''
    The lattice basis a traversal settles on, in the original lattice basis.

    If traversals 0 and i attain the same canonical form with bases C_0 and
    C_i, the linear part of the automorphism carrying one to the other is
    C_i . C_0^-1, which is unimodular because both are bases of the same
    lattice. That is how `symmetry.automorphisms` recovers the operations.

    **parameters:**
        - records: sequence
            Traversal records.

        - dim: int
            Dimension.

    **returns:**
        list or None
            `dim` integer vectors, or None when no basis exists.
    '''
    shifts = [rec[2] for rec in records]
    raw = [rec[3] for rec in records]
    _basis, transform = hermite_basis_with_transform(shifts, dim)
    if len(_basis) != dim:
        return None
    return [
        tuple(sum(row[j] * raw[j][k] for j in range(len(raw))) for k in range(dim))
        for row in transform
    ]


def _integer_edge_vectors(
    graph: PeriodicGraph,
    placement: Sequence[Sequence[Fraction]] | None = None,
) -> list[tuple[int, int, IntVec]]:
    '''
    Directed edge vectors of the barycentric drawing, cleared of denominators.

    The barycentric placement is rational, so the edge vectors are rational
    too. Multiplying every vector by the least common multiple of the
    denominators is a uniform positive scaling: it leaves every lexicographic
    comparison and every linear dependence unchanged, while allowing the whole
    traversal to run in integer arithmetic.

    **parameters:**
        - graph: PeriodicGraph
            Quotient graph, already reduced onto its own lattice.

        - placement: sequence, optional
            Barycentric placement of `graph`, when the caller has already
            solved for it. Solved for here when omitted.

    **returns:**
        tuple
            (denominator, directed edge vectors as (tail, head, integer
            vector)). The denominator is needed downstream: a cycle of the
            barycentric drawing telescopes to the integer cycle translation, so
            a record translation is `denominator * adj(B)` applied to a lattice
            vector, and the lattice the records span carries that same factor.
    '''
    placement = placement if placement is not None else barycentric_placement(graph)
    rational = edge_vectors(graph, placement)
    denominator = 1
    for (_t, _h, vec) in rational:
        for x in vec:
            denominator = denominator * x.denominator // gcd(denominator, x.denominator)
    return denominator, [
        (t, h, tuple(int(x * denominator) for x in vec))
        for (t, h, vec) in rational
    ]


def canonical_form(
    graph: PeriodicGraph,
    check_stability: bool = True,
    return_winners: bool = False,
):
    '''
    Canonical representative of the isomorphism class of `graph`.

    The graph is first reduced onto its own translation lattice and then onto
    its primitive cell, so that neither the choice of cell nor a rank-deficient
    or index-greater-than-one lattice can influence the answer. Both steps are
    needed and they are different: `PeriodicGraph.reduced` normalises the
    lattice, while `primitive.primitive` detects vertex orbits related by a
    translation, which is how an ordinary 1x1x2 supercell presents itself.

    **parameters:**
        - graph: PeriodicGraph
            Quotient graph with a connected quotient. Use
            `PeriodicGraph.components()` first if it may be disconnected.

        - check_stability: bool
            When True, refuse to canonicalise a net that is not locally
            stable rather than returning a key that depends on tie-breaking.

        - return_winners: bool
            When True, also return the reduced graph and every traversal that
            attained the canonical form. Two traversals reaching the same
            canonical representation differ by an automorphism of the net, so
            this list is the raw material of the symmetry computation
            (Delgado-Friedrichs, GD 2003, LNCS 2912, section 5: "whenever two
            d-tuples of directed edges lead to the same representation, an
            automorphism of the periodic graph has been found and all
            automorphisms must occur in this way").

    **returns:**
        PeriodicGraph
            Canonical representative; its `key_string()` is the canonical key.
            When `return_winners` is set, a triple of the canonical graph, the
            reduced graph and the list of winning traversals.

    **raises:**
        UnstableNetError
            If the net is not locally stable and `check_stability` is set.
        ValueError
            If the quotient graph is disconnected or has no valid basis.
    '''
    reduced = primitive(graph.reduced()).reduced()
    dim = reduced.dim
    if dim == 0:
        raise ValueError("the graph has no periodicity and is not a net")

    # The placement is a function of the reduced graph alone, and both the
    # stability check and the edge vectors need it, so it is solved once and
    # handed on. It is the most expensive rational arithmetic in the routine.
    placement = barycentric_placement(reduced)

    if check_stability:
        local = find_local_collisions(reduced, placement)
        if local:
            raise UnstableNetError(
                "net is not locally stable: neighbours collide in the "
                "barycentric placement, so no canonical key exists",
                collisions=local,
            )

    _denominator, directed = _integer_edge_vectors(reduced, placement)
    raw_shifts: dict[tuple[int, int, IntVec], IntVec] = {}
    for (u, v, shift) in reduced.edges:
        raw_shifts[(u, v)] = shift
    incidences: dict[int, list[tuple[int, IntVec, IntVec]]] = {
        v: [] for v in range(reduced.n_vertices)
    }
    index = 0
    for (u, v, shift) in reduced.edges:
        vec = directed[index][2]
        incidences[u].append((v, vec, shift))
        incidences[v].append((u, directed[index + 1][2], tuple(-x for x in shift)))
        index += 2

    best: tuple[Record, ...] | None = None
    finalists: list[tuple[int, list[tuple[int, int, IntVec, IntVec]], dict[int, int]]] = []
    for (start, vectors) in _candidate_bases(directed, dim):
        scaled = _scaled_inverse(vectors, dim)
        if scaled is None:
            continue
        walked = _traverse(start, scaled[0], incidences, reduced.n_vertices, dim)
        if walked is None:
            continue
        records, numbering = walked
        finished = _finalise_all(records, dim, best)
        if finished is None:
            continue
        key, verdict = finished
        if best is None or verdict < 0:
            best = key
            finalists = [(start, records, numbering)]
        elif verdict == 0:
            finalists.append((start, records, numbering))

    if best is None:
        raise ValueError("no non-singular basis of edge vectors could be found")

    winners: list[tuple[int, tuple[IntVec, ...], dict[int, int]]] = []
    if return_winners:
        for (start, records, numbering) in finalists:
            basis = finalisation_basis(records, dim)
            if basis is not None:
                winners.append((start, tuple(basis), numbering))
        if not winners:
            # Unreachable unless `hermite_basis` and its transform-tracking
            # counterpart disagree on the rank of the same generators, which
            # they are tested not to. Guarded rather than assumed, because the
            # alternative is an IndexError raised later inside `symmetry`.
            raise ValueError(
                "the canonical form was reached but no traversal yielded a "
                "lattice basis, so no automorphism frame is available"
            )

    canonical = PeriodicGraph.build(
        dim=dim,
        n_vertices=reduced.n_vertices,
        edges=[(a - 1, b - 1, shift) for (a, b, shift) in best],
        name=graph.name,
    )
    if return_winners:
        return canonical, reduced, winners
    return canonical


def canonical_key(graph: PeriodicGraph, check_stability: bool = True) -> str:
    '''
    Canonical key of a net.

    **parameters:**
        - graph: PeriodicGraph
            Quotient graph with a connected quotient.

        - check_stability: bool
            Refuse unstable nets when True.

    **returns:**
        str
            Whitespace-separated integer string.
    '''
    return canonical_form(graph, check_stability=check_stability).key_string()


__all__ = ["canonical_form", "canonical_key", "UnstableNetError"]
