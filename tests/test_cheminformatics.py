#!/usr/bin/python
'''
Tests for the choice of cheminformatics toolkit.

openbabel is the default because the identifiers shipped with the package,
the IUPAC name database among them, were generated with it. rdkit stands in
when openbabel is not installed, and what has to hold for that to be useful
is not that the two agree on everything - they do not, each has its own
canonical SMILES - but that a name looked up through one is the name looked
up through the other. That works because the database is keyed on InChIKey,
which the IUPAC algorithm defines rather than the toolkit.
'''
from __future__ import annotations

import pytest

from mofstructure import filetyper
from mofstructure import mofdeconstructor as md
from mofstructure.structure import MOFstructure

needs_openbabel = pytest.mark.skipif(
    not md.HAS_OPENBABEL, reason="openbabel is not installed")
needs_rdkit = pytest.mark.skipif(
    not md.HAS_RDKIT, reason="rdkit is not installed")

DATA = "tests/test_data"


def _first_ligand(name):
    '''First organic ligand of a test structure.'''
    return MOFstructure(filename=f"{DATA}/{name}").get_ligands()[1][0]


class TestBackendChoice:
    '''Which toolkit is picked, and what happens when one is asked for.'''

    @needs_openbabel
    def test_openbabel_is_preferred(self, monkeypatch):
        monkeypatch.delenv("MOFSTRUCTURE_CHEMINFO", raising=False)
        assert md.cheminformatics_backend() == "openbabel"

    @needs_rdkit
    def test_the_environment_can_ask_for_rdkit(self, monkeypatch):
        monkeypatch.setenv("MOFSTRUCTURE_CHEMINFO", "rdkit")
        assert md.cheminformatics_backend() == "rdkit"

    def test_an_unknown_toolkit_is_refused(self, monkeypatch):
        monkeypatch.setenv("MOFSTRUCTURE_CHEMINFO", "chemdraw")
        with pytest.raises(ImportError):
            md.cheminformatics_backend()

    def test_rdkit_is_used_when_openbabel_is_absent(self, monkeypatch):
        monkeypatch.delenv("MOFSTRUCTURE_CHEMINFO", raising=False)
        monkeypatch.setattr(md, "HAS_OPENBABEL", False)
        monkeypatch.setattr(md, "HAS_RDKIT", True)
        assert md.cheminformatics_backend() == "rdkit"

    def test_neither_toolkit_is_an_explicit_error(self, monkeypatch):
        monkeypatch.delenv("MOFSTRUCTURE_CHEMINFO", raising=False)
        monkeypatch.setattr(md, "HAS_OPENBABEL", False)
        monkeypatch.setattr(md, "HAS_RDKIT", False)
        with pytest.raises(ImportError, match="openbabel or rdkit"):
            md.cheminformatics_backend()


@needs_rdkit
class TestRdkitPerceivesFragments:
    '''
    A linker cut from its metal is an anion, and rdkit has to be told which.
    '''

    def test_terephthalate_is_the_dianion(self):
        smi = md.compute_cheminformatic_from_rdkit(
            _first_ligand("RUBTAK01.cif"))[0]
        # -2 is terephthalate; rdkit also accepts -4, which puts double bonds
        # on the ring and is not terephthalate at all.
        assert smi is not None
        assert smi.count("[O-]") == 2

    def test_trimesate_is_the_trianion(self):
        smi = md.compute_cheminformatic_from_rdkit(
            _first_ligand("HKUST-1.cif"))[0]
        assert smi is not None
        assert smi.count("[O-]") == 3

    def test_an_unperceivable_smiles_is_not_an_exception(self):
        # rdkit will not kekulize an aromatic radical, where openbabel does
        assert md._saturate_open_valences_rdkit("[n]1ccnc1") is None
        assert md.name_lookup_keys("") == []


@needs_openbabel
@needs_rdkit
class TestTheBackendsAgreeOnNames:
    '''
    The property that makes the fallback worth having.
    '''

    @pytest.mark.parametrize("cif,expected", [
        ("HKUST-1.cif", "benzene-1,3,5-triic acid"),
        ("RUBTAK01.cif", "terephthalic acid"),
    ])
    def test_the_same_ligand_gets_the_same_name(self, cif, expected,
                                                monkeypatch):
        ligand = _first_ligand(cif)
        names = filetyper.load_iupac_names()

        monkeypatch.setenv("MOFSTRUCTURE_CHEMINFO", "openbabel")
        from_openbabel = md.lookup_iupac_name(
            md.compute_openbabel_cheminformatic(ligand)[0], names)

        monkeypatch.setenv("MOFSTRUCTURE_CHEMINFO", "rdkit")
        from_rdkit = md.lookup_iupac_name(
            md.compute_cheminformatic_from_rdkit(ligand)[0], names)

        assert from_openbabel == expected
        assert from_rdkit == expected

    def test_the_inchikey_is_the_key_that_matches(self, monkeypatch):
        ligand = _first_ligand("HKUST-1.cif")
        monkeypatch.setenv("MOFSTRUCTURE_CHEMINFO", "openbabel")
        openbabel_keys = md.name_lookup_keys(
            md.compute_openbabel_cheminformatic(ligand)[0])
        monkeypatch.setenv("MOFSTRUCTURE_CHEMINFO", "rdkit")
        rdkit_keys = md.name_lookup_keys(
            md.compute_cheminformatic_from_rdkit(ligand)[0])

        assert openbabel_keys[0] == rdkit_keys[0]      # InChIKey
        assert openbabel_keys[1] != rdkit_keys[1]      # canonical SMILES
