#!/usr/bin/python
from __future__ import print_function
import pytest
from .load_test import get_test_data
import mofstructure.mofdeconstructor as MOF_deconstructor

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
    connected_components, atoms_indices_at_breaking_point, porphyrin_checker, all_regions = MOF_deconstructor.\
        secondary_building_units(ase_atom)
    metal_sbus, organic_sbus, building_unit_regions = MOF_deconstructor.\
        find_unique_building_units(
            connected_components,
            atoms_indices_at_breaking_point,
            ase_atom, porphyrin_checker,
            all_regions,
            cheminfo=True
        )
    return metal_sbus, organic_sbus, building_unit_regions, connected_components

def test_mof5(mof5):
    metal_sbus, organic_sbus, building_unit_regions, connected_components = sbu_data(mof5)
    assert len(metal_sbus) == 1
    assert len(organic_sbus) == 1
    assert len(building_unit_regions) == 2
    assert len(connected_components) == 32
    assert len(metal_sbus[0].info['point_of_extension']) == 6
    assert metal_sbus[0].info['inchikey'] == 'HCEZOJVGGICDKH-UHFFFAOYSA-N'
    assert metal_sbus[0].info['sbu_type'] == 'IRMOF_sbu'
    assert len(organic_sbus[0].info['point_of_extension']) == 2
    assert organic_sbus[0].info['inchikey'] == 'AIESRBVWAFETPR-UHFFFAOYSA-N'

def test_uio66(uio66):
    metal_sbus, organic_sbus, building_unit_regions, connected_components = sbu_data(uio66)
    # print (len(metal_sbus), len(organic_sbus), len(building_unit_regions), len(connected_components))
    # print (len(metal_sbus[0].info['point_of_extension']), metal_sbus[0].info['inchikey'], metal_sbus[0].info['sbu_type'],  len(organic_sbus[0].info['point_of_extension']), organic_sbus[0].info['inchikey'])
    assert len(metal_sbus) == 1
    assert len(organic_sbus) == 1
    assert len(building_unit_regions) == 2
    assert len(connected_components) == 28
    assert len(metal_sbus[0].info['point_of_extension']) == 12
    assert metal_sbus[0].info['inchikey'] == 'UGODGVQEWLJYCP-UHFFFAOYSA-N'
    assert metal_sbus[0].info['sbu_type'] == 'UIO66_sbu'
    assert len(organic_sbus[0].info['point_of_extension']) == 2
    assert organic_sbus[0].info['inchikey'] == 'AIESRBVWAFETPR-UHFFFAOYSA-N'


def test_dut8(dut8):
    metal_sbus, organic_sbus, building_unit_regions, connected_components = sbu_data(dut8)
    # print (len(metal_sbus), len(organic_sbus), len(building_unit_regions), len(connected_components))
    # print (len(metal_sbus[0].info['point_of_extension']), metal_sbus[0].info['inchikey'], metal_sbus[0].info['sbu_type'],  len(organic_sbus[0].info['point_of_extension']), organic_sbus[0].info['inchikey'], len(organic_sbus[1].info['point_of_extension']), organic_sbus[1].info['inchikey'])
    assert len(metal_sbus) == 1
    assert len(organic_sbus) == 2
    assert len(building_unit_regions) == 3
    assert len(connected_components) == 8
    assert len(metal_sbus[0].info['point_of_extension']) == 6
    assert metal_sbus[0].info['inchikey'] == 'ZCOYFNAPKIGPTI-UHFFFAOYSA-N'
    assert metal_sbus[0].info['sbu_type'] == 'paddlewheel'
    assert len(organic_sbus[0].info['point_of_extension']) == 2
    assert organic_sbus[0].info['inchikey'] == 'IMNIMPAHZVJRPE-UHFFFAOYSA-N'
    assert len(organic_sbus[1].info['point_of_extension']) == 2
    assert organic_sbus[1].info['inchikey'] == 'UZUHAWXNGWNYIS-UHFFFAOYSA-N'


# data = get_test_data()
# mof5, uio66, dut8 = data['MOF5'], data['UIO66'], data['DUT8']

# test_dut8(dut8)