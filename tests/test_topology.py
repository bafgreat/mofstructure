#!/usr/bin/python
from __future__ import print_function
import os
import pytest
from ase.io import read
from mofstructure import structure
from .load_test import get_test_data
from mofstructure.systre import identify_topology

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


