#!/usr/bin/env python3
'''
Geometric realisation of a net, and CGD output.

Identification is purely combinatorial, but a net still has to be *drawn*:
to look at, to hand to a downstream tool, and above all to record a topology
that carries no RCSR name. A net with no name is fully specified by its
canonical key, but a key is not something a crystallographer can open in a
viewer, so an embedding is written alongside it as a CGD file.

**why the barycentric placement is the right starting point.**
The equilibrium placement is not merely convenient, it is close to optimal for
this purpose. Delgado-Friedrichs & O'Keeffe point out that it minimises the
sum of squared edge lengths at fixed cell volume, and therefore avoids the
gratuitous entanglement that a random starting configuration produces; a
molecular-mechanics style relaxation started from a random guess spends its
time escaping self-threaded configurations instead of improving geometry.
The placement is also metric free, so it is determined by the topology alone.

**the metric.**
The barycentric placement supplies fractional coordinates but says nothing
about cell shape, because it is invariant under any affine transformation. The
cell is fixed instead by the symmetry: the invariant scalar product

    G = sum_i A_i^T A_i

built from the linear parts of the automorphism group is preserved by every
operation, so in a basis realising G the whole ideal symmetry group acts by
isometries (Delgado-Friedrichs 2005, Theorem 11). Taking the Cholesky factor
of G and scaling so that the mean edge length is one therefore produces an
embedding that realises the maximal symmetry the net can achieve, with a cell
determined by the topology and nothing else.

**two embeddings, for two different jobs.**
The barycentric placement minimises the sum of *squared* edge lengths, which
lets a few long edges pay for many short ones: real nets emerge from it with
the longest edge two or three times the shortest. That is exact, reproducible
and symmetric, and it is the right thing to record beside a canonical key.
It is a poor starting geometry for anything that has to *build* on the net,
where a linker must span every edge of one kind.

`ideal_embedding` gives the first and is the default. `refined_embedding`
gives the second, applying the penalty on edge-length variation and volume per
edge that Delgado-Friedrichs & O'Keeffe describe in section 6, and typically
brings a spread of 2.2 down to 1.0. Its result depends on the optimiser and so
is not reproducible to the last digit, which is why the two are kept apart and
why every embedding carries a `refined` flag saying which it is.

**references:**
    Delgado-Friedrichs, O. & O'Keeffe, M. (2003). Acta Cryst. A59, 351-360.
        doi:10.1107/S0108767303012017  (section 6)
    Delgado-Friedrichs, O. (2005). Discrete Comput. Geom. 33, 67-81.
        doi:10.1007/s00454-004-1147-x  (Theorem 11)
'''
from __future__ import annotations

import math
from collections.abc import Sequence

from mofstructure.graph_net.barycentric import barycentric_placement
from mofstructure.graph_net.periodic_graph import PeriodicGraph
from mofstructure.graph_net.symmetry import automorphisms, invariant_metric


def _cell_parameters(lattice: Sequence[Sequence[float]]) -> tuple[float, ...]:
    '''
    Convert lattice vectors to the six conventional cell parameters.

    **parameters:**
        - lattice: sequence
            Three lattice vectors as rows.

    **returns:**
        tuple
            (a, b, c, alpha, beta, gamma) with angles in degrees.
    '''
    def norm(v):
        return math.sqrt(sum(x * x for x in v))

    def angle(u, v):
        dot = sum(x * y for x, y in zip(u, v))
        cos = max(-1.0, min(1.0, dot / (norm(u) * norm(v))))
        return math.degrees(math.acos(cos))

    a, b, c = lattice[0], lattice[1], lattice[2]
    return (norm(a), norm(b), norm(c), angle(b, c), angle(a, c), angle(a, b))


def ideal_embedding(graph: PeriodicGraph, scale: str = "mean") -> dict[str, object]:
    '''
    Symmetric geometric realisation of a net.

    **parameters:**
        - graph: PeriodicGraph
            Net with a connected quotient graph.

        - scale: str
            `mean` normalises the mean edge length to one, `min` the shortest.

    **returns:**
        python dictionary
            Keys `lattice` (three vectors as rows), `cell` (the six
            parameters), `positions` (fractional, one per vertex orbit),
            `edges` (the quotient edges), `edge_lengths`, and
            `edge_length_spread`, the ratio of longest to shortest edge, which
            is the honest measure of how far this embedding is from one with
            uniform edges.

    **raises:**
        ValueError
            If the net is not 3-periodic, since the cell construction and the
            CGD format both assume three dimensions.
    '''
    import numpy as np

    reduced = graph.reduced()
    if reduced.dim != 3:
        raise ValueError(
            f"ideal_embedding requires a 3-periodic net, got dim={reduced.dim}"
        )

    placement = barycentric_placement(reduced)
    positions = [tuple(float(x) for x in p) for p in placement]

    group = automorphisms(reduced)
    metric = np.array(invariant_metric(group["rotations"], 3), dtype=float)
    lattice = np.linalg.cholesky(metric)

    lengths = []
    for (u, v, shift) in reduced.edges:
        frac = np.array(
            [positions[v][k] + shift[k] - positions[u][k] for k in range(3)]
        )
        lengths.append(float(np.linalg.norm(frac @ lattice)))

    if not lengths:
        raise ValueError("the net has no edges")
    reference = (sum(lengths) / len(lengths)) if scale == "mean" else min(lengths)
    lattice = lattice / reference
    lengths = [x / reference for x in lengths]

    return {
        "lattice": [[float(x) for x in row] for row in lattice],
        "cell": _cell_parameters(lattice),
        "positions": positions,
        "edges": reduced.edges,
        "edge_lengths": lengths,
        "edge_length_spread": max(lengths) / min(lengths),
        "coordination_numbers": reduced.degrees(),
        "refined": False,
    }


def refined_embedding(
    graph: PeriodicGraph,
    scale: str = "mean",
    iterations: int = 4000,
) -> dict[str, object]:
    '''
    Realisation with edges as nearly equal as the net allows.

    The barycentric placement minimises the sum of *squared* edge lengths,
    which lets a few long edges pay for many short ones: real nets come out
    of it with the longest edge two or three times the shortest. That is
    exact and reproducible and is the right thing to store, but it is a poor
    starting geometry for anything that builds on the net, where a linker has
    to span every edge of one kind. Delgado-Friedrichs & O'Keeffe deal with
    the same problem in section 6 by refining the equilibrium placement
    against a penalty on edge-length variation and on volume per edge.

    The same penalty is used here, minimised over the free coordinates and
    the lattice together. What is deliberately *not* done is to let the
    result stand in for the exact embedding: the outcome of any local
    optimiser depends on where it started and when it stopped, so two runs
    need not agree to the last digit, and `refined` is set on the result so a
    caller can tell the two apart. Store the exact one beside the key; build
    with this one.

    **parameters:**
        - graph: PeriodicGraph
            Net with a connected quotient graph.

        - scale: str
            Normalisation of the final edge lengths, as `ideal_embedding`.

        - iterations: int
            Cap on optimiser evaluations.

    **returns:**
        python dictionary
            As `ideal_embedding`, with `refined` set True and
            `edge_length_spread` reduced. `spread_before` records what the
            exact embedding gave, so the improvement is visible.

    **raises:**
        ValueError
            If the net is not 3-periodic.
    '''
    import numpy as np
    from scipy.optimize import minimize

    start = ideal_embedding(graph, scale=scale)
    positions = np.array(start["positions"], dtype=float)
    lattice = np.array(start["lattice"], dtype=float)
    edges = start["edges"]
    n_vertices = len(positions)

    def unpack(vector):
        '''Split the optimiser's vector into positions and lattice.'''
        coords = vector[: n_vertices * 3].reshape(n_vertices, 3)
        cell = vector[n_vertices * 3:].reshape(3, 3)
        return coords, cell

    def penalty(vector):
        '''
        Edge-length variation plus a term resisting collapse.

        The first term is the variance of the edge lengths, which is what
        uniformity means. The second is the total edge length divided by the
        cube root of the cell volume, without which the minimiser would
        shrink the cell to nothing and make every edge equally short.
        '''
        coords, cell = unpack(vector)
        volume = abs(np.linalg.det(cell))
        if volume < 1e-9:
            return 1e6
        lengths = np.array(
            [
                np.linalg.norm(
                    (coords[v] + np.array(shift) - coords[u]) @ cell
                )
                for (u, v, shift) in edges
            ]
        )
        if lengths.min() < 1e-9:
            return 1e6
        mean = lengths.mean()
        variation = float(((lengths - mean) ** 2).sum()) / (mean ** 2)
        density = float(lengths.sum()) / (volume ** (1.0 / 3.0))
        return variation + 0.01 * density

    # The first vertex is pinned: the placement is only defined up to a
    # translation, so leaving it free adds a flat direction the optimiser
    # would wander along without improving anything.
    guess = np.concatenate([positions.reshape(-1), lattice.reshape(-1)])
    result = minimize(
        penalty,
        guess,
        method="Nelder-Mead",
        options={"maxfev": iterations, "xatol": 1e-6, "fatol": 1e-9},
    )
    coords, cell = unpack(result.x)
    coords = coords - coords[0]

    lengths = [
        float(np.linalg.norm((coords[v] + np.array(shift) - coords[u]) @ cell))
        for (u, v, shift) in edges
    ]
    reference = (sum(lengths) / len(lengths)) if scale == "mean" else min(lengths)
    cell = cell / reference
    lengths = [x / reference for x in lengths]

    # A refinement that made things worse is not worth returning.
    if max(lengths) / min(lengths) >= start["edge_length_spread"]:
        start["spread_before"] = start["edge_length_spread"]
        return start

    return {
        "lattice": [[float(x) for x in row] for row in cell],
        "cell": _cell_parameters(cell),
        "positions": [tuple(float(x) for x in p) for p in coords],
        "edges": edges,
        "edge_lengths": lengths,
        "edge_length_spread": max(lengths) / min(lengths),
        "spread_before": start["edge_length_spread"],
        "coordination_numbers": start["coordination_numbers"],
        "refined": True,
    }


def to_cgd(
    graph: PeriodicGraph,
    name: str | None = None,
    key: str | None = None,
    rcsr: str | None = None,
    refine: bool = False,
) -> str:
    '''
    Write a net as a CGD `CRYSTAL` block.

    The block is emitted in space group P1 with every vertex orbit listed
    explicitly. That is a deliberate choice over quoting the ideal space group
    and an asymmetric unit: the orbit representatives produced here are in the
    net's own primitive setting, which need not be one of the conventional
    settings the symbol implies, and a CGD file whose coordinates and group do
    not agree is worse than one with no symmetry claimed at all. The ideal
    group is recorded in a comment instead, where it cannot be misread as a
    generator instruction.

    **parameters:**
        - graph: PeriodicGraph
            Net to write.

        - name: str or None
            Name for the NAME record.

        - key: str or None
            Canonical key, written as a comment so an unnamed net stays
            identifiable from its own file.

        - rcsr: str or None
            Archive name when the net is a named one: an RCSR symbol for a
            MOF or COF, an IZA framework-type code for a zeolite.

    **returns:**
        str
            CGD text.
    '''
    from mofstructure.graph_net.symmetry import space_group

    # Two embeddings serve two purposes. The exact one is reproducible and
    # belongs beside the key; the refined one has edges of nearly equal
    # length and is what anything building on the net wants. Which was used
    # is written into the file, because a geometry that cannot be reproduced
    # should say so rather than be mistaken for one that can.
    data = refined_embedding(graph) if refine else ideal_embedding(graph)
    label = name or graph.name or "net"
    positions = data["positions"]
    cell = data["cell"]

    try:
        group = space_group(graph)
        ideal = group.get("international")
        order = group.get("order")
    except Exception:  # noqa: BLE001 - symmetry is advisory in the output
        ideal, order = None, None

    lines = ["# Generated by mofstructure.graph_net"]
    if key:
        lines.append(f"# canonical key: {key}")
    lines.append(
        f"# net name: {rcsr}" if rcsr
        else "# net name: UNKNOWN (net not in archive)"
    )
    if ideal:
        lines.append(f"# ideal symmetry: {ideal} (order {order} modulo translations)")
    lines.append(
        "# embedding: "
        + ("refined for uniform edges, not reproducible bit for bit"
           if data.get("refined") else "exact barycentric, reproducible")
    )
    lines.append(
        "# edge length spread (max/min): "
        f"{data['edge_length_spread']:.4f}"
    )
    lines.append("CRYSTAL")
    lines.append(f"  NAME {label}")
    lines.append("  GROUP P1")
    lines.append(
        "  CELL {:.6f} {:.6f} {:.6f} {:.4f} {:.4f} {:.4f}".format(*cell)
    )
    for index, (position, degree) in enumerate(
        zip(positions, data["coordination_numbers"]), start=1
    ):
        wrapped = tuple(x - math.floor(x) for x in position)
        lines.append(
            "  NODE {} {}  {:.6f} {:.6f} {:.6f}".format(index, degree, *wrapped)
        )
    for (u, v, shift) in data["edges"]:
        tail = tuple(positions[u][k] - math.floor(positions[u][k]) for k in range(3))
        head = tuple(
            positions[v][k] + shift[k] - math.floor(positions[u][k]) for k in range(3)
        )
        lines.append(
            "  EDGE  {:.6f} {:.6f} {:.6f}   {:.6f} {:.6f} {:.6f}".format(*tail, *head)
        )
    lines.append("END")
    return "\n".join(lines) + "\n"


__all__ = ["ideal_embedding", "to_cgd", "_cell_parameters"]
