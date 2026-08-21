#!/usr/bin/env python3
'''
Local and metric-free invariants of a periodic net.

These are the descriptors that carry the weight when a net is *not* in the
RCSR archive. A canonical key answers "is this the same net as that one"; it
says nothing a chemist can read. The coordination sequence, the point symbol
and the vertex symbol are the standard published characterisation of a net,
and are what a paper reports for a topology with no three-letter symbol.

All of them are computed on the infinite covering graph, whose vertices are
pairs (v, t) of a vertex orbit and a translation, and never on the finite
quotient. Distances in the quotient are meaningless for this purpose: two
neighbours of a vertex may be the same orbit reached by different
translations.

**coordination sequence.**
For a vertex v, the k-th term n_k is the number of vertices at graph distance
exactly k from v. The sequence begins with n_1 = the coordination number.

**topological density.**
Following the Reticular Chemistry Structure Resource, the cumulative sum

    TD_n = 1 + sum_{k=1..n} n_k

counts the vertices within n shells, the vertex itself included. TD10 is the
usual quoted value and is averaged over the vertex orbits when a net has more
than one.

**cycles, rings and strong rings.**
These three are distinct and the distinction matters. A *cycle* is any closed
walk with no repeated vertex. A *ring* is a cycle that is not the sum of two
shorter cycles, equivalently a cycle on which no two vertices are joined by a
path shorter than the shortest path between them along the cycle. A *strong
ring* is a cycle that is not the sum of any number of smaller cycles. The
graph of a cube contains six 4-rings and also 6-cycles that are the sum of
three 4-rings: those 6-cycles are rings but not strong rings.
Delgado-Friedrichs & O'Keeffe recommend focusing on strong rings when
discussing the topology of a crystal graph, and note that the ring criterion
is the local one implemented here.

**point symbol and vertex symbol.**
An n-coordinated vertex has n(n-1)/2 angles, one per pair of incident edges.

  * the *point symbol*, or Schlaefli symbol, records the size of the shortest
    *cycle* through each angle, collected as A^a.B^b... with A < B < ...
  * the *vertex symbol*, or long symbol, records the size of the shortest
    *ring* through each angle, with a subscript counting how many rings of
    that size pass through the angle.

So dia has point symbol 6^6 and vertex symbol 6_2.6_2.6_2.6_2.6_2.6_2, while
pcu has point symbol 4^12.6^3 and a vertex symbol with three angles carrying
no ring at all, written with an asterisk.

**caps.**
Ring perception is unbounded in principle and a large-pore framework has
enormous rings, so `max_ring_size` bounds every search. An angle whose
shortest ring exceeds the cap is reported as `*`, exactly as the literature
writes an angle with no ring. A second cap, `path_budget`, bounds the number
of candidate paths examined at one angle; the pcu net is the standard example
of why it is needed, since its three opposite angles carry no ring at all and
a naive search enumerates every simple path up to the size cap before saying
so. Both caps are deliberate, documented approximations: without them a
single MOF pore makes the computation intractable.

**references:**
    Delgado-Friedrichs, O. & O'Keeffe, M. (2005). J. Solid State Chem. 178,
        2480-2485. doi:10.1016/j.jssc.2005.06.011  (sections 3 and 4)
    O'Keeffe, M. & Hyde, S. T. (1997). Zeolites 19, 370-374.
        (vertex symbols)
    Goetzke, R. & Klein, H.-J. (1991). J. Non-Cryst. Solids 127, 215-220.
        (rings versus strong rings)
    O'Keeffe, M., Peskov, M. A., Ramsden, S. J. & Yaghi, O. M. (2008).
        Acc. Chem. Res. 41, 1782-1789. doi:10.1021/ar800124u  (RCSR, TD10)
'''
from __future__ import annotations

from collections import Counter
from itertools import combinations
from collections.abc import Sequence

from mofstructure.graph_net.periodic_graph import PeriodicGraph

Node = tuple[int, tuple[int, ...]]

DEFAULT_MAX_RING = 12
DEFAULT_SHELLS = 10
DEFAULT_PATH_BUDGET = 20000


def _incidence_map(graph: PeriodicGraph) -> dict[int, list[tuple[int, tuple[int, ...]]]]:
    '''
    Incidences of the quotient graph, cached for repeated traversal.

    **parameters:**
        - graph: PeriodicGraph
            Net to describe.

    **returns:**
        python dictionary
            Vertex -> list of (neighbour, shift).
    '''
    return graph.incidences()


def _neighbours(
    inc: dict[int, list[tuple[int, tuple[int, ...]]]],
    node: Node,
    dim: int,
) -> list[Node]:
    '''
    Neighbours of a vertex of the infinite covering graph.

    **parameters:**
        - inc: python dictionary
            Quotient incidences.

        - node: tuple
            (vertex orbit, translation).

        - dim: int
            Dimension.

    **returns:**
        list
            Neighbouring covering-graph vertices.
    '''
    v, shift = node
    return [
        (w, tuple(shift[k] + s[k] for k in range(dim)))
        for (w, s) in inc[v]
    ]


def _bfs(
    inc: dict[int, list[tuple[int, tuple[int, ...]]]],
    source: Node,
    dim: int,
    max_depth: int,
    blocked: Node | None = None,
) -> dict[Node, int]:
    '''
    Breadth-first distances in the covering graph, optionally avoiding a vertex.

    Blocking a vertex is how a cycle through a given angle is found: the
    shortest cycle through the angle u - v - w has length 2 plus the shortest
    path from u to w that does not run back through v.

    **parameters:**
        - inc: python dictionary
            Quotient incidences.

        - source: tuple
            Starting covering-graph vertex.

        - dim: int
            Dimension.

        - max_depth: int
            Largest distance explored.

        - blocked: tuple or None
            Covering-graph vertex that may not be entered.

    **returns:**
        python dictionary
            Mapping covering-graph vertex -> distance from source.
    '''
    distance = {source: 0}
    frontier = [source]
    depth = 0
    while frontier and depth < max_depth:
        depth += 1
        nxt: list[Node] = []
        for node in frontier:
            for nb in _neighbours(inc, node, dim):
                if nb in distance or nb == blocked:
                    continue
                distance[nb] = depth
                nxt.append(nb)
        frontier = nxt
    return distance


class _DistanceCache:
    '''
    Memoised breadth-first distances in the covering graph.

    Ring perception asks for the distance between many pairs of vertices that
    lie on the same small set of cycles, so the same breadth-first search is
    requested again and again. Caching by (source, depth, blocked) turns the
    dominant cost of `vertex_symbol` from quadratic in the number of candidate
    cycles into a handful of searches.

    **parameters:**
        - inc: python dictionary
            Quotient incidences.

        - dim: int
            Dimension.
    '''

    def __init__(self, inc, dim):
        self._inc = inc
        self._dim = dim
        self._cache: dict[tuple[Node, int, Node | None], dict[Node, int]] = {}

    def distances(
        self,
        source: Node,
        max_depth: int,
        blocked: Node | None = None,
    ) -> dict[Node, int]:
        '''
        Distances from `source`, reusing a deeper cached search when possible.

        **parameters:**
            - source: tuple
                Starting covering-graph vertex.

            - max_depth: int
                Largest distance required.

            - blocked: tuple or None
                Vertex that may not be entered.

        **returns:**
            python dictionary
                Mapping covering-graph vertex -> distance.
        '''
        for depth in range(max_depth, max_depth + 4):
            hit = self._cache.get((source, depth, blocked))
            if hit is not None:
                return hit
        result = _bfs(self._inc, source, self._dim, max_depth, blocked)
        self._cache[(source, max_depth, blocked)] = result
        return result


def coordination_sequence(
    graph: PeriodicGraph,
    vertex: int = 0,
    shells: int = DEFAULT_SHELLS,
) -> list[int]:
    '''
    Coordination sequence of one vertex orbit.

    **parameters:**
        - graph: PeriodicGraph
            Net to describe.

        - vertex: int
            Vertex orbit index.

        - shells: int
            Number of shells to compute.

    **returns:**
        list
            [n_1, n_2, ..., n_shells].
    '''
    inc = _incidence_map(graph)
    source = (vertex, (0,) * graph.dim)
    distance = _bfs(inc, source, graph.dim, shells)
    counts = Counter(distance.values())
    return [counts.get(k, 0) for k in range(1, shells + 1)]


def topological_density(
    graph: PeriodicGraph,
    shells: int = DEFAULT_SHELLS,
) -> dict[str, object]:
    '''
    Cumulative topological density TD_n, per orbit and averaged.

    **parameters:**
        - graph: PeriodicGraph
            Net to describe.

        - shells: int
            Number of shells; 10 gives the conventional TD10.

    **returns:**
        python dictionary
            Keys `per_vertex`, `mean` and `shells`.
    '''
    per_vertex = []
    for v in range(graph.n_vertices):
        sequence = coordination_sequence(graph, v, shells)
        per_vertex.append(1 + sum(sequence))
    mean = sum(per_vertex) / len(per_vertex) if per_vertex else 0.0
    return {"per_vertex": per_vertex, "mean": mean, "shells": shells}


def _angles(
    inc: dict[int, list[tuple[int, tuple[int, ...]]]],
    vertex: int,
    dim: int,
) -> list[tuple[Node, Node]]:
    '''
    The n(n-1)/2 angles at a vertex, as pairs of neighbour images.

    **parameters:**
        - inc: python dictionary
            Quotient incidences.

        - vertex: int
            Vertex orbit index.

        - dim: int
            Dimension.

    **returns:**
        list
            Pairs of covering-graph neighbours of (vertex, 0).
    '''
    zero = (0,) * dim
    neighbours = [(w, s) for (w, s) in inc[vertex]]
    images = [(w, tuple(zero[k] + s[k] for k in range(dim))) for (w, s) in neighbours]
    return list(combinations(images, 2))


def _shortest_cycle_at_angle(
    inc: dict[int, list[tuple[int, tuple[int, ...]]]],
    centre: Node,
    left: Node,
    right: Node,
    dim: int,
    max_size: int,
) -> int | None:
    '''
    Size of the shortest cycle through one angle.

    **parameters:**
        - inc: python dictionary
            Quotient incidences.

        - centre: tuple
            The vertex at the apex of the angle.

        - left: tuple
            First neighbour.

        - right: tuple
            Second neighbour.

        - dim: int
            Dimension.

        - max_size: int
            Largest cycle size searched.

    **returns:**
        int or None
            Cycle size, or None when none exists within the cap.
    '''
    distance = _bfs(inc, left, dim, max_size - 2, blocked=centre)
    if right not in distance:
        return None
    return distance[right] + 2


def point_symbol(
    graph: PeriodicGraph,
    vertex: int = 0,
    max_size: int = DEFAULT_MAX_RING,
) -> str:
    '''
    Point symbol (Schlaefli symbol) of a vertex orbit.

    Records the shortest *cycle* through each of the n(n-1)/2 angles, grouped
    as A^a.B^b with increasing size. An angle carrying no cycle within the cap
    contributes an asterisk.

    **parameters:**
        - graph: PeriodicGraph
            Net to describe.

        - vertex: int
            Vertex orbit index.

        - max_size: int
            Largest cycle size searched.

    **returns:**
        str
            For example `4^12.6^3` for pcu and `6^6` for dia.
    '''
    inc = _incidence_map(graph)
    centre = (vertex, (0,) * graph.dim)
    sizes: list[int | None] = []
    for (left, right) in _angles(inc, vertex, graph.dim):
        sizes.append(
            _shortest_cycle_at_angle(inc, centre, left, right, graph.dim, max_size)
        )
    counts = Counter(s for s in sizes if s is not None)
    parts = [
        f"{size}^{count}" if count > 1 else f"{size}"
        for size, count in sorted(counts.items())
    ]
    missing = sum(1 for s in sizes if s is None)
    if missing:
        parts.append(f"*^{missing}" if missing > 1 else "*")
    return ".".join(parts)


def _is_ring(
    cache: _DistanceCache,
    cycle: Sequence[Node],
    dim: int,
) -> bool:
    '''
    Whether a cycle is a ring in the sense of Goetzke & Klein.

    A cycle is a ring when no two of its vertices are joined by a path shorter
    than the shortest path between them measured along the cycle; equivalently
    when it is not the sum of two shorter cycles.

    **parameters:**
        - cache: _DistanceCache
            Memoised distance oracle.

        - cycle: sequence
            Covering-graph vertices in cyclic order, each appearing once.

        - dim: int
            Dimension.

    **returns:**
        bool
    '''
    size = len(cycle)
    half = size // 2
    for i, node in enumerate(cycle):
        distance = cache.distances(node, half)
        for j in range(i + 1, size):
            along = min(j - i, size - (j - i))
            if along > half:
                continue
            other = cycle[j]
            actual = distance.get(other)
            if actual is not None and actual < along:
                return False
    return True


def _paths_of_length(
    inc: dict[int, list[tuple[int, tuple[int, ...]]]],
    start: Node,
    target: Node,
    length: int,
    dim: int,
    blocked: Node,
    to_target: dict[Node, int],
) -> list[list[Node]]:
    '''
    All simple paths of exactly `length` edges from `start` to `target`.

    Pruned with a precomputed lower bound on the remaining distance, so a
    branch is abandoned as soon as the target is out of reach.

    **parameters:**
        - inc: python dictionary
            Quotient incidences.

        - start: tuple
            Path start.

        - target: tuple
            Path end.

        - length: int
            Exact number of edges.

        - dim: int
            Dimension.

        - blocked: tuple
            Vertex that may not be used.

        - to_target: python dictionary
            Distances to `target`, used as the pruning bound.

    **returns:**
        list
            List of paths, each a list of covering-graph vertices.
    '''
    out: list[list[Node]] = []

    def walk(node: Node, path: list[Node], remaining: int) -> None:
        if remaining == 0:
            if node == target:
                out.append(list(path))
            return
        for nb in _neighbours(inc, node, dim):
            if nb == blocked or nb in path:
                continue
            bound = to_target.get(nb)
            if bound is None or bound > remaining - 1:
                continue
            path.append(nb)
            walk(nb, path, remaining - 1)
            path.pop()

    walk(start, [start], length)
    return out


def _rings_at_angle(
    cache: _DistanceCache,
    inc: dict[int, list[tuple[int, tuple[int, ...]]]],
    centre: Node,
    left: Node,
    right: Node,
    dim: int,
    max_size: int,
    path_budget: int = DEFAULT_PATH_BUDGET,
) -> tuple[int | None, int]:
    '''
    Size and multiplicity of the shortest rings through one angle.

    Cycles are examined in order of increasing size, and the first size at
    which a ring exists is reported together with how many rings of that size
    pass through the angle. The shortest cycle at an angle is not always a
    ring, which is precisely why the vertex symbol carries more information
    than the point symbol.

    **parameters:**
        - inc: python dictionary
            Quotient incidences.

        - centre: tuple
            Apex of the angle.

        - left: tuple
            First neighbour.

        - right: tuple
            Second neighbour.

        - dim: int
            Dimension.

        - max_size: int
            Largest ring size searched.

    **returns:**
        tuple
            (ring size or None, count).
    '''
    to_target = cache.distances(right, max_size, blocked=centre)
    shortest = to_target.get(left)
    if shortest is None:
        return (None, 0)
    spent = 0
    for path_length in range(shortest, max_size - 1):
        paths = _paths_of_length(
            inc, left, right, path_length, dim, centre, to_target
        )
        spent += len(paths)
        rings = 0
        for path in paths:
            cycle = [centre] + path
            if len(set(cycle)) != len(cycle):
                continue
            if _is_ring(cache, cycle, dim):
                rings += 1
        if rings:
            return (path_length + 2, rings)
        if spent > path_budget:
            return (None, 0)
    return (None, 0)


def vertex_symbol(
    graph: PeriodicGraph,
    vertex: int = 0,
    max_size: int = DEFAULT_MAX_RING,
) -> str:
    '''
    Vertex symbol (long symbol) of a vertex orbit.

    Records the shortest *ring* through each angle with a subscript counting
    the rings of that size, angles sorted by increasing ring size. An angle
    with no ring inside the cap is written `*`, following the convention of
    O'Keeffe & Hyde.

    Note that the special ordering convention for 4-coordinated nets, in which
    the six angles are grouped into three pairs of opposite angles, is not
    applied: which angles are opposite is a property of the coordination
    figure rather than of the graph, so it cannot be recovered from the
    topology alone. The multiset of entries is unaffected.

    **parameters:**
        - graph: PeriodicGraph
            Net to describe.

        - vertex: int
            Vertex orbit index.

        - max_size: int
            Largest ring size searched.

    **returns:**
        str
            For example `10_5.10_5.10_5` for srs.
    '''
    inc = _incidence_map(graph)
    cache = _DistanceCache(inc, graph.dim)
    centre = (vertex, (0,) * graph.dim)
    entries: list[tuple[int, str]] = []
    for (left, right) in _angles(inc, vertex, graph.dim):
        size, count = _rings_at_angle(
            cache, inc, centre, left, right, graph.dim, max_size
        )
        if size is None:
            entries.append((max_size + 1, "*"))
        elif count > 1:
            entries.append((size, f"{size}_{count}"))
        else:
            entries.append((size, f"{size}"))
    entries.sort(key=lambda item: item[0])
    return ".".join(text for (_size, text) in entries)


def describe(
    graph: PeriodicGraph,
    shells: int = DEFAULT_SHELLS,
    max_size: int = DEFAULT_MAX_RING,
) -> dict[str, object]:
    '''
    Full metric-free description of a net.

    This is what should be reported when the canonical key is absent from the
    RCSR archive: the net has no name, but it is completely characterised for
    a reader by these invariants together with its key.

    **parameters:**
        - graph: PeriodicGraph
            Net to describe.

        - shells: int
            Shells for the coordination sequence.

        - max_size: int
            Cap for cycle and ring searches.

    **returns:**
        python dictionary
            Coordination numbers, sequences, TD, point and vertex symbols,
            genus, periodicity and interpenetration multiplicity.
    '''
    return {
        "periodicity": graph.periodicity(),
        "n_vertices": graph.n_vertices,
        "n_edges": graph.n_edges,
        "genus": graph.genus(),
        "interpenetration": graph.covering_multiplicity(),
        "coordination_numbers": graph.degrees(),
        "coordination_sequences": [
            coordination_sequence(graph, v, shells) for v in range(graph.n_vertices)
        ],
        "topological_density": topological_density(graph, shells),
        "point_symbols": [
            point_symbol(graph, v, max_size) for v in range(graph.n_vertices)
        ],
        "vertex_symbols": [
            vertex_symbol(graph, v, max_size) for v in range(graph.n_vertices)
        ],
    }


__all__ = [
    "coordination_sequence",
    "topological_density",
    "point_symbol",
    "vertex_symbol",
    "describe",
]
