#!/usr/bin/python
from __future__ import print_function
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
    linker =sbu[1][0]
    assert metal_sbu.info['sbu_type'] == 'paddlewheel'
    assert metal_sbu.info['inchikey'] == 'ZCOYFNAPKIGPTI-UHFFFAOYSA-N'
    assert linker.info['inchikey'] == 'IMNIMPAHZVJRPE-UHFFFAOYSA-N'


