import csv
import json

from mofstructure.scripts import collect_data, fingerprint


SAMPLE_FINGERPRINT = {
    "clusters": {"O8Zr6": {"count": "1"}},
    "ligands": {"C8H4O4": {"count": "6"}},
    "terminal": {},
    "refinement": [],
    "fingerprint_hash": "a" * 64,
    "cluster_units": 4,
}


def test_fingerprint_cli_writes_json_and_csv(tmp_path, monkeypatch):
    cif = tmp_path / "example.cif"
    cif.write_text("test", encoding="utf-8")
    json_path = tmp_path / "fingerprints.json"
    csv_path = tmp_path / "fingerprints.csv"
    monkeypatch.setattr(
        fingerprint, "_compute_fingerprint", lambda path: SAMPLE_FINGERPRINT
    )

    result = fingerprint.main([
        str(cif), "--json", str(json_path), "--csv", str(csv_path)
    ])

    records = json.loads(json_path.read_text(encoding="utf-8"))
    record = next(iter(records.values()))
    with csv_path.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle))
    assert result == 0
    assert record["fingerprint_hash"] == "a" * 64
    assert row["fingerprint_hash"] == "a" * 64
    assert row["cluster_units"] == "4"


def test_database_collector_backfills_fingerprint(tmp_path, monkeypatch):
    cif = tmp_path / "example.cif"
    cif.write_text("test", encoding="utf-8")
    structure_data = tmp_path / "database" / "Structure_Data"
    structure_data.mkdir(parents=True)
    # An existing SBU record makes this a backfill-only run.
    (structure_data / "sbus_and_linkers.json").write_text(
        json.dumps({"example": {"n_metal_sbus": 1}}), encoding="utf-8"
    )

    class FakeMOF:
        def __init__(self, filename):
            self.filename = filename

        def get_ligand_cluster_fingerprint(self):
            return SAMPLE_FINGERPRINT

    monkeypatch.setattr(collect_data.structure, "MOFstructure", FakeMOF)
    collect_data.compile_data([str(cif)], str(tmp_path / "database"))

    records = json.loads(
        (structure_data / "fingerprint_data.json").read_text(encoding="utf-8")
    )
    assert records["example"]["fingerprint_hash"] == "a" * 64
    assert (structure_data / "fingerprint_data.csv").exists()
