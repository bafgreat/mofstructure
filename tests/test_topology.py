#!/usr/bin/python
from __future__ import print_function
import pytest
from .load_test import get_test_data
from mofstructure.systre import identify_topology

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
    assert topology == 'sql' #pcu


