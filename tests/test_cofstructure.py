#!/usr/bin/python
'''
Tests for COF deconstruction and topology.

The linkage finders are tested on small molecules built from SMILES, because a
linkage is a local pattern and a molecule is enough to show that the pattern is
matched and, just as importantly, that a look-alike is not. The topology is then
tested on two periodic structures: COF-1, which is boroxine linked, and an
idealised imine COF built here, which covers the linkage the finders are most
likely to be asked about.
'''

import os

import numpy as np
import pytest
from ase import Atoms
from ase.io import read

from mofstructure import cofstructure
from mofstructure.generate_cgd import TopologyExtractor
from mofstructure.topology import analyse

TEST_DATA = os.path.join(os.path.dirname(__file__), "test_data")


def atoms_from_smiles(smiles, pad=12.0):
    '''
    Build a 3D molecule from SMILES in a box big enough to keep it isolated.

    The box is periodic so that bond lengths can be measured with the minimum
    image convention, exactly as they are for a real structure.
    '''
    from openbabel import pybel

    mol = pybel.readstring("smi", smiles)
    mol.addh()
    mol.make3D(forcefield="mmff94", steps=250)

    positions = np.array([atom.coords for atom in mol.atoms], dtype=float)
    symbols = [pybel.ob.GetSymbol(atom.atomicnum) for atom in mol.atoms]
    positions = positions - positions.min(axis=0) + pad / 2.0
    cell = positions.max(axis=0) + pad / 2.0
    return Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)


def build_imine_hcb(a=24.08, c=20.0):
    '''
    Build an idealised 2D imine COF of 1,3,5-triformylbenzene and
    p-phenylenediamine.

    Two tritopic nodes and three ditopic linkers per hexagonal cell, which is
    the honeycomb net. The geometry is idealised rather than relaxed: bond
    lengths are set to standard values and every ring is regular, which is all
    the perception needs.
    '''
    cell = np.array([[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]])
    height = c / 2
    symbols, positions = [], []

    def add(symbol, position):
        symbols.append(symbol)
        positions.append(np.asarray(position, dtype=float))

    def unit(degrees):
        radians = np.radians(degrees)
        return np.array([np.cos(radians), np.sin(radians), 0.0])

    def add_node(centre, arms):
        for angle in arms:
            along, across = unit(angle), unit(angle + 90)
            add('C', centre + 1.40 * along)
            add('C', centre + 2.87 * along)
            add('H', centre + 2.87 * along + 1.09 * across)
        for angle in arms:
            along = unit(angle + 60)
            add('C', centre + 1.40 * along)
            add('H', centre + 2.49 * along)

    def add_linker(centre, axis):
        for sign in (1, -1):
            along = unit(axis)
            add('N', centre + sign * 2.80 * along)
            add('C', centre + sign * 1.40 * along)
        for offset in (60, 120, 240, 300):
            across = unit(axis + offset)
            add('C', centre + 1.40 * across)
            add('H', centre + 2.49 * across)

    node_a = np.array([0.0, 0.0, height])
    add_node(node_a, [90, 210, 330])
    add_node(node_a + 0.5774 * a * unit(90), [270, 30, 150])
    for axis in (90, 330, 210):
        add_linker(node_a + 0.2887 * a * unit(axis), axis)

    return Atoms(symbols=symbols, positions=np.array(positions),
                 cell=cell, pbc=True)


def linkage_counts(atoms, **options):
    '''Number of bonds found per linkage.'''
    found = cofstructure.cof_linkage_bonds(atoms, **options)
    return {name: len(bonds) for name, bonds in found.items()}


@pytest.fixture(scope="module")
def cof1():
    return read(os.path.join(TEST_DATA, "AA.cif"))


@pytest.fixture(scope="module")
def imine_cof():
    return build_imine_hcb()


# --- linkages that must be found -------------------------------------------

@pytest.mark.parametrize("name, smiles, expected", [
    ("imine", r"c1ccccc1/C=N/c1ccccc1", {"imine": 1}),
    ("ketimine", r"c1ccccc1C(=Nc1ccccc1)c1ccccc1", {"imine": 1}),
    ("hydrazone", r"c1ccccc1/C=N/NC(=O)c1ccccc1", {"imine": 1}),
    ("azine", r"c1ccccc1/C=N/N=C/c1ccccc1", {"imine": 2}),
    ("azobenzene", r"c1ccccc1/N=N/c1ccccc1", {"azo": 1}),
    ("benzanilide", r"c1ccccc1C(=O)Nc1ccccc1", {"amide": 1}),
    ("triphenyltriazine",
     r"c1ccc(cc1)-c1nc(-c2ccccc2)nc(-c2ccccc2)n1", {"triazine": 3}),
    ("boronate_ester", r"c1ccc2c(c1)OB(O2)c1ccccc1", {"boron": 1}),
    ("triphenylboroxine",
     r"c1ccc(cc1)B1OB(OB(O1)c1ccccc1)c1ccccc1", {"boron": 3}),
    ("knoevenagel", r"c1ccccc1/C=C(\C#N)c1ccccc1", {"olefin": 1}),
])
def test_linkage_is_found(name, smiles, expected):
    assert linkage_counts(atoms_from_smiles(smiles)) == expected


@pytest.mark.parametrize("name, smiles, expected", [
    # The tautomerised Tp core: three arms, so three C-N bonds to cut.
    ("tp_core",
     r"O=C1C(=CNc2ccccc2)C(=O)C(=CNc2ccccc2)C(=O)C1=CNc1ccccc1",
     {"ketoenamine": 3}),
    ("single_arm", r"O=C1CCCC(=O)C1=CNc1ccccc1", {"ketoenamine": 1}),
])
def test_ketoenamine_is_found(name, smiles, expected):
    '''
    The beta-ketoenamine of the Tp family is a C-N single bond on a
    hydrogen-bearing nitrogen, so no search for a C=N can see it.
    '''
    assert linkage_counts(atoms_from_smiles(smiles)) == expected


# --- look-alikes that must be left alone ------------------------------------

@pytest.mark.parametrize("name, smiles", [
    # An aromatic carbon next to a ring nitrogen has three neighbours, one
    # nitrogen and one hydrogen, which is the imine pattern exactly. Cutting it
    # would tear open every pyridine, bipyridine and porphyrin in a monomer.
    ("pyridine", r"c1ccncc1"),
    ("bipyridine", r"c1ccnc(c1)-c1ccccn1"),
    ("porphyrin", r"c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3"),
    # Amines are not imines.
    ("diphenylamine", r"c1ccccc1Nc1ccccc1"),
    ("triphenylamine", r"c1ccccc1N(c1ccccc1)c1ccccc1"),
    # A ketone is not an amide.
    ("benzophenone", r"c1ccccc1C(=O)c1ccccc1"),
    # An unreacted boronic acid has no ring, so there is no linkage yet.
    ("phenylboronic_acid", r"OB(O)c1ccccc1"),
])
def test_look_alike_is_not_cut(name, smiles):
    assert linkage_counts(atoms_from_smiles(smiles)) == {}


def test_pyridine_ring_survives_an_imine_cut():
    '''An imine carrying a pyridine loses the imine only.'''
    atoms = atoms_from_smiles(r"c1cc(ncc1)/C=N/c1ccccc1")
    assert linkage_counts(atoms) == {"imine": 1}


def test_hydrazide_amide_is_kept():
    '''
    The C(=O)-N of a hydrazone COF belongs inside the hydrazide monomer, so only
    the C=N is a linkage.
    '''
    found = cofstructure.cof_linkage_bonds(
        atoms_from_smiles(r"c1ccccc1/C=N/NC(=O)c1ccccc1")
    )
    assert "amide" not in found


def test_stilbene_needs_the_nitrile_rule_relaxed():
    '''
    A stilbene-cored monomer is indistinguishable from an olefin linkage, so the
    nitrile of a Knoevenagel condensation is required by default.
    '''
    atoms = atoms_from_smiles(r"c1ccccc1/C=C/c1ccccc1")
    assert linkage_counts(atoms) == {}
    assert linkage_counts(atoms, olefin_requires_nitrile=False) == {"olefin": 1}


# --- the linkage selection API ----------------------------------------------

def test_linkages_can_be_restricted():
    atoms = atoms_from_smiles(r"c1ccccc1/C=N/c1ccccc1")
    assert linkage_counts(atoms, linkages=["boron"]) == {}
    assert linkage_counts(atoms, linkages=["imine"]) == {"imine": 1}


def test_unknown_linkage_is_rejected():
    atoms = atoms_from_smiles(r"c1ccccc1")
    with pytest.raises(ValueError, match="Unknown COF linkage"):
        cofstructure.cof_linkage_bonds(atoms, linkages=["imide"])


def test_find_cn_double_bonds_shim():
    '''The pre-existing entry point still answers, now with the ring guard.'''
    assert cofstructure.find_CN_double_bonds(atoms_from_smiles(r"c1ccncc1")) == []
    assert len(cofstructure.find_CN_double_bonds(
        atoms_from_smiles(r"c1ccccc1/C=N/c1ccccc1"))) == 1


# --- ring perception --------------------------------------------------------

def test_small_rings_ignores_the_pore_ring(cof1):
    '''
    COF-1 closes a 14-membered cycle through the framework as well as its
    6-membered boroxine. Counting the large one as a ring would make every
    boron share a ring with the aryl carbon it should be cut from.
    '''
    perception = cofstructure.perceive(cof1)
    assert perception.rings
    assert all(len(ring) <= cofstructure.MAX_CHEMICAL_RING_SIZE
               for ring in perception.rings)

    boron = next(i for i, s in enumerate(perception.symbols) if s == 'B')
    aryl = next(j for j in perception.heavy[boron]
                if perception.symbols[j] == 'C')
    assert perception.in_ring(boron)
    assert not perception.share_ring(boron, aryl)


# --- deconstruction and topology --------------------------------------------

def test_cof1_deconstructs_into_boroxine_and_phenylene(cof1):
    components, bonds_to_break, _, _, breaking_pairs = \
        cofstructure.secondary_building_units(cof1)

    assert len(bonds_to_break) == 12
    assert len(breaking_pairs) == len(bonds_to_break)
    assert all(len(pair) == 5 for pair in breaking_pairs)

    formulas = sorted(cof1[c].get_chemical_formula() for c in components)
    assert formulas == ['B3O3'] * 4 + ['C6H4'] * 6


def test_cof1_topology(cof1):
    assert analyse(cof1, method="cof")['topology'] == 'hcb'


def test_cof1_through_the_extractor(cof1):
    text = TopologyExtractor(ase_atoms=cof1).build_cgd(method="cof", name="cof1")
    assert "PERIODIC_GRAPH" in text
    assert "ID cof1" in text


def test_imine_cof_deconstructs_into_node_and_linker(imine_cof):
    assert linkage_counts(imine_cof) == {"imine": 6}

    components, _, _, _, _ = cofstructure.secondary_building_units(imine_cof)
    formulas = sorted(imine_cof[c].get_chemical_formula() for c in components)
    assert formulas == ['C6H4N2'] * 3 + ['C9H6'] * 2


def test_imine_cof_topology(imine_cof):
    '''
    The tritopic nodes stay vertices and the ditopic linkers are spliced into
    edges, which is what makes the net hcb rather than its subdivision.
    '''
    edges, node_atoms = cofstructure.cof_topology_graph(imine_cof)
    assert len(node_atoms) == 2
    assert len(edges) == 3
    assert analyse(imine_cof, method="cof")['topology'] == 'hcb'


def test_ditopic_linkers_are_kept_when_asked(imine_cof):
    edges, node_atoms = cofstructure.cof_topology_graph(
        imine_cof, collapse_ditopic=False
    )
    assert len(node_atoms) == 5
    assert len(edges) == 6


def test_structure_without_a_supported_linkage_says_so(cof1):
    with pytest.raises(ValueError, match="no supported COF linkage"):
        cofstructure.cof_topology_graph(cof1, linkages=["imine"])
