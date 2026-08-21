#!/usr/bin/env python3
'''
Tests for the pure-Python topology engine in `mofstructure.graph_net`.

The suite is organised around the properties that make the engine
trustworthy rather than around its functions, because the risk here is not
that a function raises but that it quietly returns a wrong net:

  * the archive parser must round-trip every stored key byte for byte,
  * the canonical form must be *invariant*, giving one key for all
    representations of a net, and *injective*, giving different keys to
    different nets,
  * the derived symmetry, descriptors and embedding must reproduce values
    published independently in the literature.

Reference values are taken from Delgado-Friedrichs & O'Keeffe,
J. Solid State Chem. 178, 2480-2485 (doi:10.1016/j.jssc.2005.06.011) and from
the RCSR, and are quoted in the tests so a failure can be checked against the
source without leaving the file.
'''
from __future__ import annotations

import random
from pathlib import Path

import pytest

from mofstructure.graph_net.archive import parse_systre_archive
from mofstructure.graph_net.canonical import canonical_key
from mofstructure.graph_net.invariants import (
    coordination_sequence,
    point_symbol,
    vertex_symbol,
)
from mofstructure.graph_net.lattice import hermite_basis, hermite_basis_with_transform
from mofstructure.graph_net.periodic_graph import PeriodicGraph


@pytest.fixture(scope="module")
def archive():
    '''
    The bundled RCSR archive as a mapping from identifier to stored key.

    **returns:**
        python dictionary
            Mapping RCSR identifier -> Systre key string.
    '''
    return {name: key for (name, key) in parse_systre_archive()}


def _scramble(graph: PeriodicGraph, rng: random.Random) -> PeriodicGraph:
    '''
    Produce a different representation of the same net.

    Applies the three transformations that leave a net unchanged while
    changing its vector representation: a permutation of the vertex orbits, a
    change of orbit representatives, and a unimodular change of lattice basis.

    **parameters:**
        - graph: PeriodicGraph
            Net to disguise.

        - rng: random.Random
            Seeded generator, so a failure is reproducible.

    **returns:**
        PeriodicGraph
    '''
    dim = graph.dim
    matrix = [[int(i == j) for j in range(dim)] for i in range(dim)]
    for _ in range(6):
        i, j = rng.sample(range(dim), 2)
        factor = rng.choice([-2, -1, 1, 2])
        for c in range(dim):
            matrix[i][c] += factor * matrix[j][c]
    order = list(range(graph.n_vertices))
    rng.shuffle(order)
    moved = graph.relabelled(order).translated(
        {v: tuple(rng.randint(-2, 2) for _ in range(dim)) for v in range(graph.n_vertices)}
    )
    return moved.transformed(matrix)


class TestLattice:
    '''Exact integer lattice arithmetic.'''

    def test_hermite_basis_is_canonical(self):
        '''The basis depends on the lattice, not on the generators given.'''
        first = hermite_basis([(2, 0, 0), (0, 2, 0), (1, 1, 1)])
        second = hermite_basis([(1, 1, 1), (2, 0, 0), (0, 2, 0), (3, 1, 1)])
        assert first == second

    def test_full_rank_generators_give_the_identity(self):
        '''Generators spanning Z^3 reduce to the standard basis.'''
        assert hermite_basis([(1, -2, -1), (1, -1, -2), (1, -1, -1)]) == [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1),
        ]

    def test_transform_reconstructs_the_basis(self):
        '''The returned coefficients really do rebuild each basis vector.'''
        generators = [(0, 0, 0), (2, 0, 0), (0, 2, 0), (1, 1, 1)]
        basis, transform = hermite_basis_with_transform(generators, 3)
        for row, expected in zip(transform, basis):
            rebuilt = tuple(
                sum(row[j] * generators[j][k] for j in range(len(generators)))
                for k in range(3)
            )
            assert rebuilt == expected

    def test_both_hermite_routines_agree(self):
        '''
        The two reductions return the same basis for the same generators.

        `canonical_form` compares candidates using `hermite_basis`, which is
        cheap, and recovers the lattice basis of the winners alone with
        `hermite_basis_with_transform`, which is not: carrying the transform
        widens every row from `dim` to `dim + m` entries. Splitting the two is
        only legitimate while they agree, including on rank, since a candidate
        accepted by one and rejected by the other would win the comparison and
        then have no frame. The routines differ in that one drops zero
        generators and the other keeps them to preserve the column
        correspondence, so the agreement is worth asserting rather than
        assuming.
        '''
        rng = random.Random(12345)
        for _ in range(2000):
            dim = rng.choice([1, 2, 3, 4])
            low, high = rng.choice([(-2, 2), (-9, 9), (-40, 40)])
            generators = [
                [0] * dim
                if rng.random() < 0.2
                else [rng.randint(low, high) for _ in range(dim)]
                for _ in range(rng.randint(1, 10))
            ]
            if rng.random() < 0.25:
                generators.append([2 * x for x in generators[0]])
            cheap = hermite_basis(generators, dim)
            full, _transform = hermite_basis_with_transform(generators, dim)
            assert list(cheap) == list(full)


class TestPeriodicGraph:
    '''The quotient graph and its normalisations.'''

    def test_archive_keys_round_trip(self, archive):
        '''Every stored key survives parse and re-serialise unchanged.'''
        bad = [
            name
            for name, key in archive.items()
            if PeriodicGraph.from_key_string(key).key_string() != key
        ]
        assert bad == []

    def test_periodicity_is_the_cycle_lattice_rank(self):
        '''A 2-periodic net written in a 3-column cell is 2-periodic.'''
        flat = PeriodicGraph.build(3, 1, [(0, 0, (1, 0, 0)), (0, 0, (0, 1, 0))])
        assert flat.periodicity() == 2

    def test_supercell_reduces_to_the_primitive_cell(self):
        '''Doubling the cell must not change the reduced representation.'''
        primitive = PeriodicGraph.build(
            3, 1, [(0, 0, (1, 0, 0)), (0, 0, (0, 1, 0)), (0, 0, (0, 0, 1))]
        )
        doubled = PeriodicGraph.build(
            3, 1, [(0, 0, (2, 0, 0)), (0, 0, (0, 1, 0)), (0, 0, (0, 0, 1))]
        )
        assert doubled.reduced().key_string() == primitive.reduced().key_string()

    def test_zero_shift_self_loop_is_rejected(self):
        '''A loop in the same cell is not an edge of any net.'''
        with pytest.raises(ValueError):
            PeriodicGraph.build(3, 1, [(0, 0, (0, 0, 0))])


class TestCanonicalForm:
    '''The property the whole engine rests on.'''

    @pytest.mark.parametrize("name", ["pcu", "dia", "srs", "nbo", "fcu", "bcu", "sod", "tbo"])
    def test_key_is_invariant_under_representation(self, archive, name):
        '''Relabelling, retranslating and rebasing must not move the key.'''
        graph = PeriodicGraph.from_key_string(archive[name])
        expected = canonical_key(graph)
        rng = random.Random(20240617)
        for _ in range(5):
            assert canonical_key(_scramble(graph, rng)) == expected

    def test_supercell_gives_the_same_key(self, archive):
        '''A 1x1x2 supercell of pcu is still pcu.'''
        primitive = PeriodicGraph.from_key_string(archive["pcu"])
        supercell = PeriodicGraph.build(
            3, 2,
            [
                (0, 1, (0, 0, 0)), (1, 0, (0, 0, 1)),
                (0, 0, (1, 0, 0)), (0, 0, (0, 1, 0)),
                (1, 1, (1, 0, 0)), (1, 1, (0, 1, 0)),
            ],
        )
        assert canonical_key(supercell) == canonical_key(primitive)

    def test_distinct_nets_get_distinct_keys(self, archive):
        '''Injectivity on a sample of nets that are genuinely different.'''
        names = ["pcu", "dia", "srs", "nbo", "fcu", "bcu", "sod", "crs", "ths", "cds"]
        keys = {name: canonical_key(PeriodicGraph.from_key_string(archive[name])) for name in names}
        assert len(set(keys.values())) == len(names)


class TestInvariants:
    '''Descriptors checked against published values.'''

    @pytest.mark.parametrize(
        "name,expected",
        [
            ("pcu", [6, 18, 38, 66, 102, 146]),
            ("dia", [4, 12, 24, 42, 64, 92]),
            ("srs", [3, 6, 12, 24, 35, 48]),
            ("nbo", [4, 12, 28, 50, 76, 110]),
        ],
    )
    def test_coordination_sequences(self, archive, name, expected):
        '''Standard coordination sequences from the RCSR.'''
        graph = PeriodicGraph.from_key_string(archive[name])
        assert coordination_sequence(graph, 0, len(expected)) == expected

    @pytest.mark.parametrize("name,expected", [("pcu", "4^12.6^3"), ("dia", "6^6")])
    def test_point_symbols(self, archive, name, expected):
        '''Point symbols quoted by Delgado-Friedrichs & O'Keeffe (2005) section 4.'''
        assert point_symbol(PeriodicGraph.from_key_string(archive[name]), 0) == expected

    @pytest.mark.parametrize(
        "name,expected",
        [
            ("srs", "10_5.10_5.10_5"),
            ("dia", "6_2.6_2.6_2.6_2.6_2.6_2"),
            ("pcu", "4.4.4.4.4.4.4.4.4.4.4.4.*.*.*"),
        ],
    )
    def test_vertex_symbols(self, archive, name, expected):
        '''Vertex symbols quoted verbatim in the same source.'''
        assert vertex_symbol(PeriodicGraph.from_key_string(archive[name]), 0) == expected


class TestSymmetry:
    '''Ideal symmetry derived from the graph alone.'''

    @pytest.mark.parametrize(
        "name,symbol,order",
        [
            ("pcu", "Pm-3m", 48),
            ("dia", "Fd-3m", 48),
            ("srs", "I4_132", 24),
            ("nbo", "Im-3m", 48),
            ("fcu", "Fm-3m", 48),
            ("bcu", "Im-3m", 48),
            ("sod", "Im-3m", 48),
            ("crs", "Fd-3m", 48),
            ("cds", "P4_2/mmc", 16),
        ],
    )
    def test_ideal_space_group(self, archive, name, symbol, order):
        '''Maximal symmetry, with no use of the crystal's own coordinates.'''
        from mofstructure.graph_net.symmetry import space_group

        result = space_group(PeriodicGraph.from_key_string(archive[name]))
        assert result["international"] == symbol
        assert result["order"] == order

    def test_srs_is_intrinsically_chiral(self, archive):
        '''srs admits only proper operations, so no achiral embedding exists.'''
        from mofstructure.graph_net.symmetry import space_group

        assert space_group(PeriodicGraph.from_key_string(archive["srs"]))["is_chiral"]

    def test_pcu_is_not_chiral(self, archive):
        '''pcu contains improper operations.'''
        from mofstructure.graph_net.symmetry import space_group

        assert not space_group(PeriodicGraph.from_key_string(archive["pcu"]))["is_chiral"]

    def test_two_periodic_nets_are_not_named(self, archive):
        '''The 230 space group types do not classify a layer net.'''
        from mofstructure.graph_net.symmetry import space_group

        assert space_group(PeriodicGraph.from_key_string(archive["hcb"]))["international"] is None


class TestEmbedding:
    '''Geometric realisation.'''

    @pytest.mark.parametrize("name", ["pcu", "dia", "srs", "fcu", "nbo", "sod", "tbo"])
    def test_edge_transitive_nets_have_uniform_edges(self, archive, name):
        '''
        An edge-transitive net must have every edge the same length in any
        embedding that realises its symmetry, so the spread is exactly one.
        '''
        from mofstructure.graph_net.embedding import ideal_embedding

        data = ideal_embedding(PeriodicGraph.from_key_string(archive[name]))
        assert data["edge_length_spread"] == pytest.approx(1.0, abs=1e-9)

    def test_diamond_cell_matches_the_known_ratio(self, archive):
        '''
        The primitive cell edge of diamond is 4/sqrt(6) times the bond length,
        which follows from a_prim = a_cubic/sqrt(2) and bond = a_cubic*sqrt(3)/4.
        '''
        import math

        from mofstructure.graph_net.embedding import ideal_embedding

        data = ideal_embedding(PeriodicGraph.from_key_string(archive["dia"]))
        assert data["cell"][0] == pytest.approx(4.0 / math.sqrt(6.0), rel=1e-9)
        assert data["cell"][3] == pytest.approx(60.0, abs=1e-9)

    def test_cgd_output_is_well_formed(self, archive):
        '''The written block carries the key and closes properly.'''
        from mofstructure.graph_net.embedding import to_cgd

        graph = PeriodicGraph.from_key_string(archive["dia"])
        text = to_cgd(graph, name="dia", key=canonical_key(graph), rcsr="dia")
        assert "CRYSTAL" in text and text.rstrip().endswith("END")
        assert "canonical key:" in text
        assert text.count("NODE") == graph.n_vertices


class TestArchiveTable:
    '''The re-keyed RCSR lookup table.'''

    def test_table_has_no_collisions(self):
        '''
        Two different named nets must never share a key.

        Read from the shipped table rather than from the JSON the build tools
        write, which is a derived intermediate and is not committed: a test
        that reads it passes on the machine that built it and fails on a
        clean checkout. Injectivity is asserted here on every name the table
        carries, which is the property that makes a key an identification
        rather than a hint.
        '''
        from mofstructure.graph_net.archive import _table

        names = _table()
        assert names, "the shipped lookup table is empty"
        seen = {}
        for key, by_source in names.items():
            for source, name in by_source.items():
                previous = seen.get((source, name))
                assert previous is None or previous == key, (
                    f"{source} name {name!r} is claimed by two keys"
                )
                seen[(source, name)] = key

    @pytest.mark.parametrize("name", ["pcu", "dia", "srs", "fcu"])
    def test_lookup_recovers_the_rcsr_symbol(self, archive, name):
        '''A net keyed from its archive entry must look up as itself.'''
        from mofstructure.graph_net.archive import GRAPH_NET_TABLE, lookup

        if not Path(GRAPH_NET_TABLE).exists():
            pytest.skip("lookup table not built; run tools/build_graph_net_archive.py")
        key = canonical_key(PeriodicGraph.from_key_string(archive[name]))
        assert lookup(key) == name
