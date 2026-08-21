from pathlib import Path

import math
import numpy as np
from ase.io import read

from analysis.controlled_perturbation_benchmark import (
    atomic_graph_hash,
    fingerprint_payload,
    missing_hydrogen,
    missing_linker,
    reordered_atoms,
    shifted_origin,
    supercell_2x,
    wilson_interval,
)


TEST_DATA = Path(__file__).parent / "test_data"


def test_invariant_transforms_and_finite_graph_baseline():
    atoms = read(TEST_DATA / "Cr.cif")
    rng = np.random.default_rng(7)
    pristine = fingerprint_payload(atoms)["fingerprint_hash"]
    graph = atomic_graph_hash(atoms)

    for transform in (reordered_atoms, shifted_origin):
        variant = transform(atoms, rng).atoms
        assert fingerprint_payload(variant)["fingerprint_hash"] == pristine
        assert atomic_graph_hash(variant) == graph

    repeated = supercell_2x(atoms, rng).atoms
    assert fingerprint_payload(repeated)["fingerprint_hash"] == pristine
    assert atomic_graph_hash(repeated) != graph


def test_defect_transforms_change_fingerprint():
    atoms = read(TEST_DATA / "RUBTAK01.cif")
    rng = np.random.default_rng(11)
    pristine = fingerprint_payload(atoms)["fingerprint_hash"]

    linker_defect = missing_linker(atoms, rng)
    assert linker_defect.status == "created"
    assert fingerprint_payload(linker_defect.atoms)["fingerprint_hash"] != pristine

    hydrogen_defect = missing_hydrogen(atoms, rng)
    assert hydrogen_defect.status == "created"
    assert fingerprint_payload(hydrogen_defect.atoms)["fingerprint_hash"] != pristine


def test_wilson_interval_is_bounded():
    low, high = wilson_interval(50, 100)
    assert 0.39 < low < 0.41
    assert 0.59 < high < 0.61
    empty_low, empty_high = wilson_interval(0, 0)
    assert math.isnan(empty_low)
    assert math.isnan(empty_high)
