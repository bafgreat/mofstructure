#!/usr/bin/env python3
'''
Maximal symmetry of a periodic net, derived from its topology alone.

The symmetry reported here is the *ideal* symmetry: the largest
crystallographic group any embedding of the net can achieve. It is computed
from the graph, never from the coordinates of the crystal the net came from,
so a badly distorted experimental structure and its idealised counterpart give
the same answer. Delgado-Friedrichs & O'Keeffe note that this is exactly why
structures are so often published in the wrong symmetry, and that knowing the
maximal embeddable group gives a useful upper bound.

**where the operations come from.**
Every isomorphism of periodic graphs induces an affine map between their
barycentric drawings (Delgado-Friedrichs 2004, Theorem 1), and for a stable
net the automorphism group is isomorphic to a crystallographic space group
(2004, Corollary 1). The automorphisms themselves fall out of the canonical
form for nothing: two traversals of Algorithm 5 that reach the same canonical
representation differ precisely by an automorphism, and all automorphisms
arise this way. Each winning traversal carries a frame, the ordered tuple of
edge vectors it used as a coordinate basis, so if the reference traversal used
the frame B_0 and another winner used B_i, the linear part of the
corresponding automorphism is

    A_i = B_i . B_0^-1,                                                  (1)

and its translation part is fixed by where the start vertex goes,

    b_i = p(start_i) - A_i . p(start_0).                                 (2)

This is the construction of Delgado-Friedrichs & O'Keeffe (2003, section 4):
a triplet of directed edges in linearly independent directions determines an
affine mapping completely.

**why a metric is needed.**
The A_i are integer matrices in the lattice basis, but they are only
*isometries* with respect to a suitable metric. Following the standard trick
in mathematical crystallography, the invariant scalar product

    G = sum_i A_i^T A_i                                                  (3)

is positive definite and is preserved by every A_i, so in a basis realising G
the whole group acts by isometries (Delgado-Friedrichs 2005, Theorem 11). The
Cholesky factor of G supplies the lattice used for naming the group.

**naming.**
Systre hands its computed operations to `sgtbx` to recognise the group among
the 230 standard settings. The same job is done here by spglib's
`get_spacegroup_type_from_symmetry`, which takes rotations and translations
rather than atomic positions and so consumes the combinatorially derived
operations directly.

**intrinsic chirality.**
If every operation of the maximal group is proper, no embedding of the net can
be achiral, and the net is *intrinsically chiral*. The converse does not hold:
a topologically chiral embedding may still admit an achiral one, as
Delgado-Friedrichs & O'Keeffe show with their P4(1)22 and Ama2 realisations of
one CdSO4-derived net.

**what is reported when the net is only locally stable.**
`canonical_form` admits a net that is locally stable, which is weaker than the
stability Corollary 1 asks for: distinct vertices may still share a
barycentric position as long as no two neighbours of one vertex do. For such a
net the two groups genuinely differ. An automorphism that carries a vertex
onto one it collides with moves no point of the drawing, so it is realised by
the identity, and the operations found here are deduplicated on the isometry
they induce. What comes back is therefore the maximal crystallographic
symmetry, which is what an embedding can actually achieve, and `order` is the
order of that group rather than of the automorphism group, which is strictly
larger. Delgado-Friedrichs & O'Keeffe make exactly this point in section 6,
and the ladder nets of Delgado-Friedrichs, O'Keeffe & Treacy are the
systematic family where it bites: every vertex collides with a
symmetry-equivalent one, so the automorphism group is not isomorphic to any
space group at all. The two groups coincide precisely when the net is stable.

**references:**
    Delgado-Friedrichs, O. & O'Keeffe, M. (2003). Acta Cryst. A59, 351-360.
        doi:10.1107/S0108767303012017  (sections 4 and 6)
    Delgado-Friedrichs, O. (2004). Graph Drawing (GD 2003), LNCS 2912, 178-189.
        (Theorems 1 and 2, Corollary 1; section 5: automorphisms from
        repeated canonical representations)
    Delgado-Friedrichs, O. (2005). Discrete Comput. Geom. 33, 67-81.
        doi:10.1007/s00454-004-1147-x  (Theorem 11, the invariant metric)
    Delgado-Friedrichs, O., O'Keeffe, M. & Treacy, M. M. J. (2020).
        Acta Cryst. A76, 735-738. doi:10.1107/S2053273320012905  (ladder nets)
'''
from __future__ import annotations

from fractions import Fraction
from collections.abc import Sequence

from mofstructure.graph_net.barycentric import barycentric_placement
from mofstructure.graph_net.canonical import canonical_form
from mofstructure.graph_net.periodic_graph import PeriodicGraph

RatMatrix = list[list[Fraction]]


def _invert(columns: Sequence[Sequence[int]], dim: int) -> RatMatrix | None:
    '''
    Inverse of the matrix whose columns are `columns`, over the rationals.

    **parameters:**
        - columns: sequence
            `dim` integer vectors taken as columns.

        - dim: int
            Side length.

    **returns:**
        list or None
            Rows of the inverse, or None when singular.
    '''
    mat = [[Fraction(columns[c][r]) for c in range(dim)] for r in range(dim)]
    inv = [[Fraction(int(r == c)) for c in range(dim)] for r in range(dim)]
    for col in range(dim):
        pivot = next((r for r in range(col, dim) if mat[r][col] != 0), None)
        if pivot is None:
            return None
        mat[col], mat[pivot] = mat[pivot], mat[col]
        inv[col], inv[pivot] = inv[pivot], inv[col]
        scale = Fraction(1) / mat[col][col]
        mat[col] = [x * scale for x in mat[col]]
        inv[col] = [x * scale for x in inv[col]]
        for row in range(dim):
            if row == col or mat[row][col] == 0:
                continue
            factor = mat[row][col]
            mat[row] = [mat[row][k] - factor * mat[col][k] for k in range(dim)]
            inv[row] = [inv[row][k] - factor * inv[col][k] for k in range(dim)]
    return inv


def _matmul(left: Sequence[Sequence[Fraction]], right: Sequence[Sequence[Fraction]], dim: int) -> RatMatrix:
    '''
    Product of two square rational matrices.

    **parameters:**
        - left: sequence
            Rows of the left factor.

        - right: sequence
            Rows of the right factor.

        - dim: int
            Side length.

    **returns:**
        list
            Rows of the product.
    '''
    return [
        [
            sum((Fraction(left[r][k]) * Fraction(right[k][c]) for k in range(dim)), Fraction(0))
            for c in range(dim)
        ]
        for r in range(dim)
    ]


def _apply(matrix: Sequence[Sequence[Fraction]], vec: Sequence[Fraction], dim: int) -> tuple[Fraction, ...]:
    '''
    Apply a rational matrix to a rational vector.

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
    return tuple(
        sum((Fraction(matrix[r][c]) * Fraction(vec[c]) for c in range(dim)), Fraction(0))
        for r in range(dim)
    )


def automorphisms(graph: PeriodicGraph) -> dict[str, object]:
    '''
    Affine operations of the maximal symmetry group, modulo lattice translations.

    Applies equations (1) and (2) to every traversal that attained the
    canonical form.

    **parameters:**
        - graph: PeriodicGraph
            Net with a connected quotient graph.

    **returns:**
        python dictionary
            Keys `dim`, `order`, `rotations` (integer matrices in the lattice
            basis), `translations` (fractional), and `permutations` (the
            induced maps on vertex orbits).

    **raises:**
        ValueError
            If the linear parts fail to come out integral, which would mean
            the net is not stable enough for the construction to apply.
    '''
    _canonical, reduced, winners = canonical_form(graph, return_winners=True)
    dim = reduced.dim
    placement = barycentric_placement(reduced)

    start0, frame0, numbering0 = winners[0]
    inverse0 = _invert(frame0, dim)
    if inverse0 is None:
        raise ValueError("the reference frame is singular")
    inverse_by_number = {number: v for v, number in numbering0.items()}

    rotations: list[list[list[int]]] = []
    translations: list[tuple[Fraction, ...]] = []
    permutations: list[tuple[int, ...]] = []
    seen = set()

    for (start, frame, numbering) in winners:
        frame_matrix = [[Fraction(frame[c][r]) for c in range(dim)] for r in range(dim)]
        linear = _matmul(frame_matrix, inverse0, dim)
        if any(x.denominator != 1 for row in linear for x in row):
            raise ValueError(
                "an automorphism has a non-integral linear part; the net is "
                "not stable in the sense required for symmetry analysis"
            )
        shift = tuple(
            placement[start][k] - _apply(linear, placement[start0], dim)[k]
            for k in range(dim)
        )
        rotation = tuple(tuple(int(x) for x in row) for row in linear)
        translation = tuple(x - (x.numerator // x.denominator) for x in shift)
        signature = (rotation, translation)
        if signature in seen:
            continue
        seen.add(signature)
        rotations.append([list(row) for row in rotation])
        translations.append(translation)
        permutations.append(
            tuple(
                inverse_by_number[numbering[v]] if numbering[v] in inverse_by_number else v
                for v in range(reduced.n_vertices)
            )
        )

    return {
        "dim": dim,
        "order": len(rotations),
        "rotations": rotations,
        "translations": translations,
        "permutations": permutations,
    }


def invariant_metric(rotations: Sequence[Sequence[Sequence[int]]], dim: int) -> list[list[float]]:
    '''
    Positive definite metric preserved by every operation, equation (3).

    **parameters:**
        - rotations: sequence
            Integer linear parts of the group.

        - dim: int
            Dimension.

    **returns:**
        list
            The metric tensor as a nested list of floats.
    '''
    metric = [[0.0] * dim for _ in range(dim)]
    for rot in rotations:
        for r in range(dim):
            for c in range(dim):
                metric[r][c] += sum(rot[k][r] * rot[k][c] for k in range(dim))
    return metric


def space_group(graph: PeriodicGraph, symprec: float = 1e-5) -> dict[str, object]:
    '''
    Maximal crystallographic symmetry of a net, named by spglib.

    **parameters:**
        - graph: PeriodicGraph
            Net with a connected quotient graph.

        - symprec: float
            Tolerance handed to spglib. The operations are exact integers, so
            this only guards the floating point lattice built from the metric.

    **returns:**
        python dictionary
            Keys `order`, `international`, `number`, `hall`, `crystal_system`,
            `is_chiral` and `rotations`. Naming keys are None for a net that
            is not 3-periodic, since the 230 space group types only classify
            three dimensions.
    '''
    import numpy as np
    import spglib

    group = automorphisms(graph)
    dim = group["dim"]
    rotations = group["rotations"]
    translations = group["translations"]
    proper = all(_determinant(rot, dim) == 1 for rot in rotations)

    result: dict[str, object] = {
        "order": group["order"],
        "is_chiral": proper,
        "rotations": rotations,
        "international": None,
        "number": None,
        "hall": None,
        "crystal_system": None,
    }
    if dim != 3:
        return result

    metric = np.array(invariant_metric(rotations, dim), dtype=float)
    lattice = np.linalg.cholesky(metric)
    # spglib is deprecating the error handling that returns None, in favour
    # of raising. Accepting both leaves the naming keys as None either way,
    # rather than tying this to one spglib release.
    try:
        dataset = spglib.get_spacegroup_type_from_symmetry(
            np.array(rotations, dtype="intc"),
            np.array([[float(x) for x in t] for t in translations],
                     dtype="double"),
            lattice=lattice,
            symprec=symprec,
        )
    except Exception:
        dataset = None
    if dataset is not None:
        as_dict = dataset if isinstance(dataset, dict) else dataset.__dict__
        result["international"] = as_dict.get("international_short") or as_dict.get("international")
        result["number"] = as_dict.get("number")
        result["hall"] = as_dict.get("hall_symbol") or as_dict.get("hall")
        result["crystal_system"] = as_dict.get("crystal_system")
    return result


def _determinant(matrix: Sequence[Sequence[int]], dim: int) -> int:
    '''
    Determinant of a small integer matrix.

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
        minor = [[matrix[r][k] for k in range(dim) if k != c] for r in range(1, dim)]
        total += ((-1) ** c) * matrix[0][c] * _determinant(minor, dim - 1)
    return total


__all__ = ["automorphisms", "invariant_metric", "space_group"]
