from copy import deepcopy

from analysis import fingerprint_usefulness as analysis


def _record(name, count="1", refinement=None):
    record = {
        "cluster_units": 1,
        "clusters": {
            "O8Zr6": {
                "count": "1",
                "contacts": {"10": "1"},
                "capped": {"0": "1"},
            }
        },
        "ligands": {
            "C8H4O4": {
                "count": count,
                "contacts": {"2": count},
                "denticity": {"1": count},
            }
        },
        "terminal": {},
        "refinement": refinement or [["colour-a", "1"]],
        "mof_name": name,
        "source": f"/{name}.cif",
    }
    record["fingerprint_hash"] = analysis.fingerprint_hash_from_record(record)
    return record


def test_hash_verification_restores_integer_histogram_keys():
    record = _record("ABCDEF.MOF_subset")
    # The production hash sorts integer key 2 before 10. A JSON-loaded mapping
    # sorts string key "10" before "2" unless the keys are restored.
    record["clusters"]["O8Zr6"]["contacts"]["2"] = "1"
    record["fingerprint_hash"] = analysis.fingerprint_hash_from_record(record)

    assert analysis.fingerprint_hash_from_record(record) == record["fingerprint_hash"]
    assert analysis.naive_json_roundtrip_hash(record) != record["fingerprint_hash"]


def test_nested_resolution_levels_add_information():
    first = _record("ABCDEF.MOF_subset")
    duplicate = deepcopy(first)
    duplicate["mof_name"] = "ABCDEF01.MOF_subset"
    second_stoichiometry = _record("GHIJKL.MOF_subset", count="2")
    second_refinement = deepcopy(second_stoichiometry)
    second_refinement["mof_name"] = "MNOPQR.MOF_subset"
    second_refinement["refinement"] = [["colour-b", "1"]]
    second_refinement["fingerprint_hash"] = analysis.fingerprint_hash_from_record(
        second_refinement
    )
    records = {
        "a": first,
        "b": duplicate,
        "c": second_stoichiometry,
        "d": second_refinement,
    }

    bundle = analysis.analyze_records(records, top_n=5)
    resolution = {
        row["level"]: row["unique_signatures"]
        for row in bundle["resolution_rows"]
    }

    assert resolution == {
        "chemistry_presence": 1,
        "stoichiometry": 2,
        "local_coordination": 2,
        "full_refinement": 3,
    }
    assert bundle["summary"]["full_fingerprint"]["maximum_group_size"] == 2
    assert bundle["summary"]["hash_integrity"][
        "mismatches_after_restoring_integer_histogram_keys"
    ] == 0
