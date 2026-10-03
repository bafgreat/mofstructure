"""
The canonical key must not depend on how a crystal is written down.

Each test changes only the representation of one structure: the origin of the
cell, the cell itself or the order of the atoms. A key that moves under any of
these is not an identifier, and two files of one material would be counted as
two materials.

LORLEJ and EVAKAQ are a layered and a rod MOF whose all-node nets contain
edges built inside the metal SBU. An origin shift moves some of those bonds
across the cell boundary, which exposed a sign error in the translation given
to those edges. EMIYUV deconstructs to a one-periodic net whose barycentric
placement has collisions, for which neither the primitive cell nor the key is
well defined. COF_586 (CoRE COF 586) is a single hydrazone layer in a cell
small enough that one linkage joins a building unit to its own neighbouring
image, and its two-connected units are spliced into edges; it exposed a cut
that was discarded because both ends lay in one component and a sign error
in the splice.

The sbus and ligand_cluster reductions are not tested here. They collapse a
rod SBU to a single vertex, which is not a well-defined periodic graph, so
for rod MOFs their nets depend on the cell by construction.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pytest
from ase.build import make_supercell
from ase.io import read

from mofstructure.graph_net.canonical import UnstableNetError, canonical_key
from mofstructure.graph_net.periodic_graph import PeriodicGraph
from mofstructure.structure import MOFstructure
from mofstructure.topology import _identify_into, analyse

DATA = Path(__file__).resolve().parent / "test_data"


def _read(name):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return read(DATA / name)


def _key(atoms, method="all_node"):
    record = analyse(MOFstructure(ase_atoms=atoms).remove_guest(),
                     method=method, timeout=300)
    return record["status"], record.get("key")


def _variants(atoms, seed=0):
    rng = np.random.default_rng(seed)
    for k in range(3):
        moved = atoms.copy()
        moved.translate(rng.random(3) @ atoms.cell)
        moved.wrap()
        yield f"origin shift {k}", moved
    yield "permutation", atoms[rng.permutation(len(atoms))]
    yield "2x1x1 supercell", atoms.repeat((2, 1, 1))
    yield "cell setting", make_supercell(atoms, [[1, 1, 0], [0, 1, 0], [0, 0, 1]],
                                         wrap=True)


@pytest.mark.parametrize("filename", ["LORLEJ.cif", "EVAKAQ.cif", "Cr.cif"])
def test_all_node_key_is_independent_of_the_representation(filename):
    atoms = _read(filename)
    status, reference = _key(atoms)
    assert status == "ok" and reference
    for label, variant in _variants(atoms):
        assert _key(variant) == (status, reference), f"key moved under {label}"


def test_cof_key_is_independent_of_the_representation():
    atoms = _read("COF_586.cif")
    status, reference = _key(atoms, method="cof")
    assert status == "ok" and reference
    for label, variant in _variants(atoms, seed=3):
        assert _key(variant, method="cof") == (status, reference), \
            f"key moved under {label}"


def test_an_unstable_net_is_refused_in_every_cell():
    atoms = _read("EMIYUV.cif")
    statuses = {_key(variant)[0] for _, variant in _variants(atoms)}
    statuses.add(_key(atoms)[0])
    assert statuses == {"unidentified"}


def _graph_supercell(graph, axis, n=2):
    edges = []
    for copy in range(n):
        for (u, v, shift) in graph.edges:
            shift = list(shift)
            total = copy + shift[axis]
            shift[axis] = total // n
            edges.append((u + graph.n_vertices * copy,
                          v + graph.n_vertices * (total % n), tuple(shift)))
    return PeriodicGraph.build(dim=graph.dim, n_vertices=graph.n_vertices * n,
                               edges=edges)


def test_a_net_with_collisions_has_no_key():
    # Quotient graph of the EMIYUV chain. Vertices 0 and 3, and 1 and 2, share
    # a position in the barycentric placement although no vertex has two
    # neighbours at the same point, so the net is locally stable but not
    # stable. Before the global test it received a key that changed when the
    # same chain was written in a doubled cell.
    chain = PeriodicGraph.build(dim=3, n_vertices=4, edges=[
        (0, 2, (-1, 0, 0)), (0, 2, (0, 0, 0)), (0, 3, (0, 1, 0)),
        (1, 2, (0, -1, 0)), (1, 3, (0, 0, 0)), (1, 3, (1, 0, 0)),
    ])
    for graph in (chain, _graph_supercell(chain, 0), _graph_supercell(chain, 0, 3)):
        with pytest.raises(UnstableNetError):
            canonical_key(graph)


def test_a_stable_net_has_the_same_key_in_a_supercell():
    # pcu written in its primitive cell and as 2x1x1, 1x2x1 and 1x1x3 cells.
    pcu = PeriodicGraph.build(dim=3, n_vertices=1, edges=[
        (0, 0, (1, 0, 0)), (0, 0, (0, 1, 0)), (0, 0, (0, 0, 1)),
    ])
    reference = canonical_key(pcu)
    for axis, n in ((0, 2), (1, 2), (2, 3)):
        assert canonical_key(_graph_supercell(pcu, axis, n)) == reference


def test_a_timeout_is_reported_as_a_timeout():
    # MOZ takes several seconds to key, so a one-second budget expires inside
    # the canonical form. The timeout used to be caught by the per-component
    # error handling and reported as "unidentified", which reads as a property
    # of the net rather than of the budget.
    record = analyse(str(DATA / "MOZ.cif"), method="zeol", timeout=1)
    assert record["status"] == "timeout"
