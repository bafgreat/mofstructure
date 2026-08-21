from copy import deepcopy

import pytest

from analysis import fingerprint_similarity as similarity


def _record(metal="Zn", ligand="C8H4O4", contacts="2", refinement="a"):
    return {
        "cluster_units": 1,
        "clusters": {
            metal: {
                "count": "1",
                "contacts": {"4": "1"},
                "capped": {"0": "1"},
            }
        },
        "ligands": {
            ligand: {
                "count": "2",
                "contacts": {contacts: "2"},
                "denticity": {"1": "4"},
            }
        },
        "terminal": {},
        "refinement": [[refinement, "1"]],
        "fingerprint_hash": f"hash-{metal}-{ligand}-{contacts}-{refinement}",
        "mof_name": f"{metal}-{ligand}-{contacts}",
        "source": "/test.cif",
    }


def test_identical_detailed_records_score_one():
    record = _record()
    scores = similarity.detailed_similarity(
        record, deepcopy(record), similarity.WEIGHT_PRESETS["balanced"]
    )

    assert scores == {
        "composition": pytest.approx(1.0),
        "local": pytest.approx(1.0),
        "shape": pytest.approx(1.0),
        "refinement": pytest.approx(1.0),
        "overall": pytest.approx(1.0),
    }


def test_connectivity_channel_recognizes_metal_substitution():
    zinc = _record(metal="Zn", refinement="zn")
    copper_same_shape = _record(metal="Cu", refinement="cu")
    copper_different_shape = _record(metal="Cu", contacts="3", refinement="other")
    weights = similarity.WEIGHT_PRESETS["connectivity"]

    same_shape = similarity.detailed_similarity(zinc, copper_same_shape, weights)
    different_shape = similarity.detailed_similarity(
        zinc, copper_different_shape, weights
    )

    assert same_shape["shape"] == pytest.approx(1.0)
    assert different_shape["shape"] < same_shape["shape"]
    assert same_shape["overall"] > different_shape["overall"]


def test_cluster_coordination_mismatch_is_not_overwhelmed_by_ligand_counts():
    six_connected = _record()
    eighteen_connected = _record(metal="Zr", refinement="zr")
    eighteen_connected["clusters"]["Zr"]["contacts"] = {"18": "1"}
    eighteen_connected["ligands"]["C8H4O4"]["count"] = "9"
    eighteen_connected["ligands"]["C8H4O4"]["contacts"] = {"2": "9"}
    eighteen_connected["ligands"]["C8H4O4"]["denticity"] = {"1": "18"}

    scores = similarity.detailed_similarity(
        six_connected,
        eighteen_connected,
        similarity.WEIGHT_PRESETS["balanced"],
    )

    assert scores["shape"] < 0.85
    assert scores["overall"] < 0.80


def test_sparse_index_candidates_are_reranked_with_exact_scores():
    records = {
        "query": _record(),
        "exact": deepcopy(_record()),
        "analogue": _record(metal="Cu", refinement="cu"),
        "different": _record(metal="Cu", contacts="3", refinement="other"),
    }
    records["query"]["mof_name"] = "query"
    records["exact"]["mof_name"] = "exact"
    records["analogue"]["mof_name"] = "analogue"
    records["different"]["mof_name"] = "different"
    record_ids = list(records)
    weights = similarity.WEIGHT_PRESETS["balanced"]
    matrix = similarity.build_sparse_matrix(records, record_ids, 256, weights)

    results = similarity.nearest_neighbors(
        0,
        matrix,
        record_ids,
        records,
        weights,
        top=2,
        candidate_pool=4,
    )

    assert results[0]["neighbor_record_id"] == "exact"
    assert results[0]["overall_similarity"] == pytest.approx(1.0)
    assert results[1]["neighbor_record_id"] == "analogue"
