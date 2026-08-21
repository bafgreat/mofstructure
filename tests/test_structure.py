#!/usr/bin/python
from mofstructure import structure
from .load_test import get_test_data


def test_structure():
    '''
    Test to ensure that zeo++ works efficiently in computing
    geometric properties of MOFs
    '''
    data = get_test_data()
    ase_atom = data['DUT8']

    # filename = './test_data/SARSUC.cif'

    mof = structure.MOFstructure(ase_atom)

    assert len(mof.get_oms()) == 8
    cluster_ligand = mof.get_ligands()
    assert len(cluster_ligand) == 2
    assert cluster_ligand[1][0].info['inchikey'] == 'DESKXQISJOCIDX-UHFFFAOYSA-N'

    sbu = mof.get_sbu()
    metal_sbu = sbu[0][0]
    linker = sbu[1][0]
    assert metal_sbu.info['sbu_type'] == 'paddlewheel'
    assert metal_sbu.info['inchikey'] == 'ZCOYFNAPKIGPTI-UHFFFAOYSA-N'
    assert linker.info['inchikey'] == 'IMNIMPAHZVJRPE-UHFFFAOYSA-N'
    topology = mof.get_topology(method="all_node")
    assert topology.get('topology') == 'pcu'
    assert topology.get('dimension') == 3
    assert topology.get('td10') == 1561
    # The hash changed when identification moved off Systre, and deliberately.
    # It used to digest Systre's relaxed geometry, rounded coordinates and
    # cell; it now digests the canonical key. The new one identifies the net
    # rather than one drawing of it, so the same framework in a supercell or
    # with its atoms reordered hashes the same, which the old one did not.
    # The version prefix is there so a value stored under the old scheme is
    # recognisable rather than merely wrong.
    assert topology.get('key_version') == 'graph_net/1'
    assert topology.get('topology_hash') == (
        'graph_net/1:sha256:'
        'c09a19983feb32ef9c5337506781ebaab4c89debba4d1971e53489440bd1351c'
    )
    assert topology.get('topology_hash') == topology.get('key_hash')
    assert topology.get('key')
    assert topology.get('cgd')
