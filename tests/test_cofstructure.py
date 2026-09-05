#!/usr/bin/python
'''
Tests for COF linkage detection, deconstruction and topology.

Molecular fixtures cover linkage patterns and exclusions. Periodic fixtures
cover boroxine-linked COF-1, an idealised imine COF and azodioxy-linked NPN-1.
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


# --- bond perception --------------------------------------------------------

def line_of_atoms(symbols, gaps):
    '''Place atoms on a line with successive gaps in an isolated periodic box.'''
    positions, x = [], 0.0
    for gap in [0.0] + list(gaps):
        x += gap
        positions.append((x, 0.0, 0.0))
    return Atoms(symbols=symbols, positions=np.array(positions),
                 cell=np.eye(3) * 30.0, pbc=True)


def test_a_hydrogen_bond_is_not_a_covalent_bond():
    '''Remove the 1.50 A H...O contact while retaining N-H and O-C.'''
    atoms = line_of_atoms(['N', 'H', 'O', 'C'], [1.01, 1.50, 1.22])
    perception = cofstructure.perceive(atoms)
    assert perception.graph[1] == [0]
    assert perception.n_hydrogen[2] == 0


def test_a_contact_between_two_disorder_images_is_not_a_bond():
    '''Remove a 1.90 A C...C contact between two bonded pairs.'''
    atoms = line_of_atoms(['C', 'C', 'C', 'C'], [1.42, 1.90, 1.42])
    perception = cofstructure.perceive(atoms)
    assert 2 not in perception.graph[1]
    assert perception.graph[0] == [1] and perception.graph[3] == [2]


def test_a_long_bond_is_kept_when_it_is_the_only_one():
    '''Retain the long methoxy C-O bond found in TpOMe-PaNO2.'''
    atoms = line_of_atoms(['C', 'O', 'C'], [1.37, 1.67])
    perception = cofstructure.perceive(atoms)
    assert perception.heavy[2] == [1]


def test_two_hydrogens_do_not_pair_into_a_molecule():
    '''Keep each hydrogen on carbon even when the H...H contact is shorter.'''
    atoms = line_of_atoms(['C', 'H', 'H', 'C'], [1.20, 0.90, 1.20])
    perception = cofstructure.perceive(atoms)
    assert perception.graph[1] == [0]
    assert perception.graph[2] == [3]


def test_an_imine_drawn_without_hydrogens_is_still_found():
    '''Detect the same imine bond after removing explicit hydrogens.'''
    atoms = atoms_from_smiles(r"c1ccccc1/C=N/c1ccccc1")
    heavy = atoms[[a.index for a in atoms if a.symbol != 'H']]
    assert linkage_counts(heavy) == {"imine": 1}


# --- bridge fallback --------------------------------------------------------

@pytest.mark.parametrize("name, smiles, expected", [
    # Two cuts leave a ditopic C=O or C=S fragment.
    ("urea", r"c1ccccc1NC(=O)Nc1ccccc1", {"bridge": 2}),
    ("thiourea", r"c1ccccc1NC(=S)Nc1ccccc1", {"bridge": 2}),
    # This enaminone falls outside the Tp ketoenamine pattern.
    ("enaminone", r"O=C(c1ccccc1)/C=C/Nc1ccccc1", {"bridge": 1}),
    ("sulfonamide", r"c1ccccc1S(=O)(=O)Nc1ccccc1", {"bridge": 1}),
])
def test_an_unnamed_linkage_is_still_cut(name, smiles, expected):
    assert linkage_counts(atoms_from_smiles(smiles)) == expected


@pytest.mark.parametrize("name, smiles", [
    # Exclude C-C spacers within monomers.
    ("tetraphenylethylene",
     r"C(=C(c1ccccc1)c1ccccc1)(c1ccccc1)c1ccccc1"),
    ("tolane", r"c1ccccc1C#Cc1ccccc1"),
    ("biphenyl", r"c1ccccc1-c1ccccc1"),
    # A single-atom bridge has no internal bond to cut.
    ("diphenyl_ether", r"c1ccccc1Oc1ccccc1"),
    # Explicit CH2 groups exclude these glycol and polyether spacers.
    ("glycol_bridge", r"c1ccccc1OCCOc1ccccc1"),
    ("polyether_bridge", r"c1ccccc1OCCOCCOc1ccccc1"),
])
def test_the_fallback_leaves_a_monomer_alone(name, smiles):
    assert linkage_counts(atoms_from_smiles(smiles)) == {}


def test_the_fallback_defers_to_a_named_finder():
    '''Do not add bridge cuts to linkages already identified by named finders.'''
    for smiles in (r"c1ccccc1/C=N/c1ccccc1",
                   r"c1ccccc1C(=O)Nc1ccccc1",
                   r"c1ccccc1/N=N/c1ccccc1"):
        assert "bridge" not in cofstructure.cof_linkage_bonds(
            atoms_from_smiles(smiles)
        )


def test_a_bridge_needs_two_different_cores():
    '''Treat directly bonded phenyl rings as one core.'''
    atoms = atoms_from_smiles(r"c1ccccc1-c1ccccc1")
    perception = cofstructure.perceive(atoms)
    assert cofstructure.acyclic_bridges(perception) == []
    assert len(set(cofstructure.ring_systems(perception).values())) == 1


@pytest.mark.parametrize("name, smiles", [
    # Terminal oxygen substituents must not prevent the N-N cut.
    ("azoxybenzene", r"c1ccccc1/N=[N+](\[O-])c1ccccc1"),
    ("azodioxy", r"c1ccccc1[N+]([O-])=[N+]([O-])c1ccccc1"),
])
def test_an_oxidised_azo_is_still_an_azo(name, smiles):
    assert linkage_counts(atoms_from_smiles(smiles)) == {"azo": 1}


def test_a_nitro_group_is_not_an_azo():
    '''Exclude nitro nitrogen, which has no nitrogen neighbour.'''
    assert linkage_counts(atoms_from_smiles(r"c1ccccc1[N+](=O)[O-]")) == {}


@pytest.mark.parametrize('shift_periodic_image', [False, True])
def test_npn1_retains_oxygen_on_its_building_units_and_has_dia_topology(
    shift_periodic_image,
):
    '''The periodic azodioxy framework must cut N-N while retaining N-O.'''
    atoms = read(os.path.join(TEST_DATA, "NPN-1.cif"))
    if shift_periodic_image:
        nitrogen = next(atom.index for atom in atoms if atom.symbol == 'N')
        atoms.positions[nitrogen] += atoms.cell[0]
    assert linkage_counts(atoms) == {"azo": 4}
    components, cuts, _, _, breaking_pairs = \
        cofstructure.secondary_building_units(atoms)
    assert all(atoms[i].symbol == atoms[j].symbol == 'N' for i, j in cuts)
    assert sorted(atoms[c].get_chemical_formula() for c in components) == \
        ['C25H16N4O4', 'C25H16N4O4']
    if shift_periodic_image:
        assert any(any(pair[2:]) for pair in breaking_pairs)
    result = analyse(atoms)
    assert result['status'] == 'ok'
    assert result['material'] == 'cof'
    assert result['topology'] == 'dia'


def test_bridge_fallback_runs_after_named_finders_in_any_requested_order():
    atoms = atoms_from_smiles(r"c1ccccc1/C=N/c1ccccc1")
    assert linkage_counts(atoms, linkages=['bridge', 'imine']) == {'imine': 1}
