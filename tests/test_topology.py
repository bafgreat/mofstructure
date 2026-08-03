#!/usr/bin/python
from __future__ import print_function
import os
import numpy as np
import pytest
from ase.io import read
from mofstructure import structure
from .load_test import get_test_data
from mofstructure.systre import identify_topology
from mofstructure.generate_cgd import ligand_cluster_graph, net_geometry

TEST_DATA = os.path.join(os.path.dirname(__file__), "test_data")

@pytest.fixture(scope="module")
def data():
    return get_test_data()


@pytest.fixture(scope="module")
def mof5(data):
    return data['MOF5']


@pytest.fixture(scope="module")
def uio66(data):
    return data['UIO66']


@pytest.fixture(scope="module")
def dut8(data):
    return data['DUT8']


def sbu_data(ase_atom):
    '''
    Function to compile secondary building units and regions of MOFs.

    Parameters
    ----------
    ase_atom : ASE atoms object
    '''
    res = identify_topology(ase_atom)
    return res.topology


def test_mof5(mof5):
    topology = sbu_data(mof5)
    assert topology == 'pcu'

def test_uio66(uio66):
    topology = sbu_data(uio66)

    assert topology == 'fcu'

def test_dut8(dut8):
    topology = sbu_data(dut8)
    assert topology == 'pcu' #pcu


def _cif(name):
    return read(os.path.join(TEST_DATA, name))


def test_hkust1_tbo():
    # tritopic BTC: the polytopic linker must be its own node, giving tbo (not reo)
    mof = structure.MOFstructure(_cif("HKUST-1.cif"))
    assert mof.get_topology(method="all_node")["topology"] == "tbo"


def test_rod_all_node_vs_sbus():
    # MIL-53 rod: all_node keeps the chain (rna), sbus collapses it (pcu).
    # Matches CrystalNets AllNodes=rna.
    mof = structure.MOFstructure(_cif("Cr.cif"))
    assert mof.get_topology(method="all_node")["topology"] == "rna"
    assert mof.get_topology(method="sbus")["topology"] == "pcu"


def test_rod_single_node():
    # MIL-53 rod, single_node: metals stay separate, each linker merges to one
    # vertex. Matches CrystalNets SingleNodes=bpq.
    mof = structure.MOFstructure(_cif("Cr.cif"))
    assert mof.get_topology(method="single_node")["topology"] == "bpq"


def test_single_node_discrete_unchanged():
    # single_node must not change discrete-SBU frameworks
    assert structure.MOFstructure(_cif("HKUST-1.cif")).get_topology(
        method="single_node")["topology"] == "tbo"


def test_ligand_cluster_hkust1_incidence_net():
    atoms = _cif("HKUST-1.cif")
    edges, ligand_nodes, node_atoms = ligand_cluster_graph(atoms)
    metal_nodes = set(node_atoms) - ligand_nodes

    assert len(metal_nodes) == 24
    assert len(ligand_nodes) == 32
    assert all(
        (u in metal_nodes and v in ligand_nodes)
        or (v in metal_nodes and u in ligand_nodes)
        for u, v, *_ in edges
    )

    degree = {node: 0 for node in node_atoms}
    for u, v, *_ in edges:
        degree[u] += 1
        degree[v] += 1
    assert {degree[node] for node in metal_nodes} == {4}
    assert {degree[node] for node in ligand_nodes} == {3}
    assert all(
        sum(atoms[i].symbol == "C" for i in node_atoms[node]) == 9
        and sum(atoms[i].symbol == "O" for i in node_atoms[node]) == 6
        for node in ligand_nodes
    )


def test_ligand_cluster_keeps_distinct_periodic_incidences():
    edges, ligand_nodes, _ = ligand_cluster_graph(_cif("Cr.cif"))
    incidences = [edge for edge in edges if edge[0] != edge[1]]
    assert len(incidences) == len(set(incidences))
    assert all((u in ligand_nodes) != (v in ligand_nodes) for u, v, *_ in incidences)


def test_ligand_cluster_partitions_the_atoms():
    # every atom belongs to exactly one vertex: no atom is counted both in a
    # ligand and in the cluster it coordinates
    atoms = _cif("HKUST-1.cif")
    _, ligand_nodes, node_atoms = ligand_cluster_graph(atoms)
    ligand_atoms = set().union(*(node_atoms[n] for n in ligand_nodes))
    cluster_atoms = set().union(
        *(node_atoms[n] for n in set(node_atoms) - ligand_nodes)
    )
    assert not ligand_atoms & cluster_atoms
    assert ligand_atoms | cluster_atoms == set(range(len(atoms)))


def test_ligand_cluster_rod_keeps_equivalent_linkers_equivalent():
    # MIL-53's four BDC linkers are one symmetry orbit (Imma), so contracting
    # the rod must not hand them different coordination numbers
    edges, ligand_nodes, node_atoms = ligand_cluster_graph(_cif("Cr.cif"))
    degree = {node: 0 for node in node_atoms}
    for u, v, *_ in edges:
        degree[u] += 1
        degree[v] += 1
    assert {degree[node] for node in ligand_nodes} == {2}


def test_ligand_cluster_survives_a_supercell():
    # the same crystal in a bigger box is the same framework
    from mofstructure.generate_cgd import ligand_cluster_fingerprint
    atoms = _cif("Cr.cif")
    assert (
        ligand_cluster_fingerprint(atoms)["fingerprint_hash"]
        == ligand_cluster_fingerprint(atoms.repeat((1, 1, 2)))["fingerprint_hash"]
    )


def test_ligand_cluster_ignores_coordinated_solvent():
    # a methanol on an open Cu site is not a linker: the net stays tbo and the
    # framework part of the fingerprint is untouched, but the solvent is recorded
    from ase import Atom
    from mofstructure.generate_cgd import ligand_cluster_fingerprint

    atoms = _cif("HKUST-1.cif")
    cu = next(i for i, s in enumerate(atoms.get_chemical_symbols()) if s == "Cu")
    distances = atoms.get_distances(cu, range(len(atoms)), mic=True)
    axial = -sum(
        (lambda v: v / np.linalg.norm(v))(
            atoms.get_distance(cu, j, mic=True, vector=True)
        )
        for j in np.argsort(distances)[1:6] if atoms[j].symbol == "O"
    )
    axial /= np.linalg.norm(axial)

    solvated = atoms.copy()
    solvated.append(Atom("O", atoms.positions[cu] + 2.2 * axial))
    solvated.append(Atom("C", atoms.positions[cu] + 3.63 * axial))

    assert structure.MOFstructure(solvated).get_topology(
        method="ligand_cluster")["topology"] == "tbo"

    # counts are quoted per metal cluster, so one methanol among 24 is 1/24
    fingerprint = ligand_cluster_fingerprint(solvated)
    assert fingerprint["terminal"] == {
        "CO": {"count": "1/24", "denticity": {1: "1/24"}}
    }
    assert (
        fingerprint["ligands"] == ligand_cluster_fingerprint(atoms)["ligands"]
    )


def test_ligand_cluster_fingerprint_sees_a_missing_linker():
    # a missing-linker defect leaves two clusters one contact short
    from mofstructure.generate_cgd import ligand_cluster_fingerprint
    from mofstructure import mofdeconstructor

    atoms = _cif("RUBTAK01.cif")
    pristine = ligand_cluster_fingerprint(atoms)
    assert pristine["clusters"] == {
        "O8Zr6": {"count": "1", "contacts": {12: "1"}, "capped": {0: "1"}}
    }
    assert pristine["ligands"]["C8H4O4"]["contacts"] == {2: "6"}

    components, *_ = mofdeconstructor.ligands_and_metal_clusters(atoms)
    linker = next(
        c for c in components
        if not any(atoms[int(i)].symbol == "Zr" for i in c)
    )
    defective = atoms[[
        i for i in range(len(atoms)) if i not in {int(x) for x in linker}
    ]]

    # half the clusters are now one linker short
    fingerprint = ligand_cluster_fingerprint(defective)
    assert fingerprint["clusters"]["O8Zr6"]["contacts"] == {11: "1/2", 12: "1/2"}
    assert fingerprint["fingerprint_hash"] != pristine["fingerprint_hash"]


def test_draw_topology():
    # the net drawn over the real structure has the expected node counts:
    # MIL-53 all_node = 4 metal + 8 carboxyl = 12; single_node merges to 8
    plotly = pytest.importorskip("plotly")  # noqa: F841
    mof = structure.MOFstructure(_cif("Cr.cif"))
    fig = mof.draw_topology(method="all_node", show_topology=True)
    markers = sum(
        len(t.x) for t in fig.data
        if t.mode == "markers" and t.name in ("metal nodes", "organic nodes")
    )
    # The 12 nodes in the requested cell plus any neighbouring periodic nodes
    # needed to terminate the displayed crossing edges.
    assert markers >= 12
    fig2 = mof.draw_topology(method="single_node", show_topology=True)
    assert sum(
        len(t.x) for t in fig2.data
        if t.mode == "markers" and t.name in ("metal nodes", "organic nodes")
    ) >= 8

    # It is displayed as a clean structure/network, not an axis-based plot.
    assert fig.layout.scene.xaxis.visible is False
    assert fig.layout.scene.yaxis.visible is False
    assert fig.layout.scene.zaxis.visible is False
    assert fig.layout.scene.camera.projection.type == "orthographic"
    assert fig.layout.scene.bgcolor == "#f4f7f5"
    assert fig.layout.paper_bgcolor == "#f4f7f5"
    assert fig.layout.legend.bgcolor == "rgba(255,255,255,0.82)"


def test_draw_topology_marks_both_ends_of_every_edge():
    plotly = pytest.importorskip("plotly")  # noqa: F841
    mof = structure.MOFstructure(_cif("Cr.cif"))
    fig = mof.draw_topology(
        method="all_node", show_structure=False, show_unit_cell=False,
        show_linker_sbu=False, show_topology=True,
    )
    edge_trace = next(t for t in fig.data if t.name == "topology edges")
    node_traces = [
        t for t in fig.data if t.name in ("metal nodes", "organic nodes")
    ]
    node_coords = {
        tuple(round(float(value), 7) for value in xyz)
        for trace in node_traces for xyz in zip(trace.x, trace.y, trace.z)
    }
    edge_points = [
        tuple(round(float(value), 7) for value in xyz)
        for xyz in zip(edge_trace.x, edge_trace.y, edge_trace.z)
        if xyz[0] is not None
    ]
    assert edge_points
    assert all(point in node_coords for point in edge_points)


def test_draw_topology_shows_linker_sbu_mapping_for_every_method():
    plotly = pytest.importorskip("plotly")  # noqa: F841
    mof = structure.MOFstructure(_cif("HKUST-1.cif"))
    for method in ("sbus", "all_node", "single_node", "ligand_cluster"):
        fig = mof.draw_topology(
            method=method, show_structure=False, show_unit_cell=False
        )
        names = {trace.name for trace in fig.data}
        assert "topology edges" not in names
        assert "metal nodes" not in names
        assert "organic nodes" not in names
        assert f"{method} connections" in names
        assert "SBU / metal centres" in names
        assert "organic / linker centres" in names
        contacts = next(
            trace for trace in fig.data
            if trace.name == f"{method} connections"
        )
        assert contacts.line.color == "#008000"
        assert contacts.line.width == 7
        assert contacts.line.dash == "solid"
        assert sum(
            x is None
            for trace in fig.data if trace.name == f"{method} connections"
            for x in trace.x
        ) > 0
        centre_coords = {
            tuple(round(float(value), 7) for value in xyz)
            for trace in fig.data
            if trace.name in (
                "SBU / metal centres", "organic / linker centres"
            )
            for xyz in zip(trace.x, trace.y, trace.z)
        }
        contact_endpoints = {
            tuple(round(float(value), 7) for value in xyz)
            for xyz in zip(contacts.x, contacts.y, contacts.z)
            if xyz[0] is not None
        }
        assert contact_endpoints <= centre_coords


def test_draw_geometry_preserves_periodic_edges():
    positions, _, edges, cell = net_geometry(
        _cif("Cr.cif"), method="ligand_cluster"
    )
    self_edges = [edge for edge in edges if edge[0] == edge[1]]
    assert self_edges
    assert all(tuple(edge[2:]) != (0, 0, 0) for edge in self_edges)
    assert all(
        np.linalg.norm(
            positions[v] + np.asarray((sx, sy, sz)) @ cell - positions[u]
        ) > 1e-6
        for u, v, sx, sy, sz in edges
    )


def test_single_node_drawing_has_no_zero_self_edges():
    _, _, edges, _ = net_geometry(_cif("Cr.cif"), method="single_node")
    assert all(not (u == v and (sx, sy, sz) == (0, 0, 0))
               for u, v, sx, sy, sz in edges)
