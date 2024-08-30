#!/usr/bin/python
from __future__ import print_function
import os
import pytest
from ase.io import read

test_dir = os.path.dirname(os.path.abspath(__file__))

list_of_files = [
    os.path.join(test_dir, 'test_data/AA.cif'),
    os.path.join(test_dir, 'test_data/EDUSIF.cif'),
    os.path.join(test_dir, 'test_data/RUBTAK01.cif'),
    os.path.join(test_dir, 'test_data/SARSUC.cif')
]

def load_data(list_of_test_data):
    all_data = []
    for mole_file in list_of_test_data:
        ase_atom = read(mole_file)
        all_data.append(ase_atom)
    return all_data

def get_test_data():
    ase_data = load_data(list_of_files)
    data = {'COF1': ase_data[0], 'MOF5': ase_data[1], 'UIO66': ase_data[2], 'DUT8': ase_data[3]}
    return data

@pytest.fixture
def test_data():
    return get_test_data()

def test_structure_count(test_data):
    assert len(test_data) == 4
