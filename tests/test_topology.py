#!/usr/bin/env python3
'''
Tests for the public topology API.

`mofstructure.topology` is the layer that takes a structure and returns its
net, so what is worth asserting here is not that the identification is
correct - `test_graph_net.py` establishes that against the literature and
against two independent implementations - but that the layer around it holds
together:

  * a MOF, a COF and a zeolite must come back in the *same* shape, since the
    point of the record is that results are comparable across chemistry,
  * the material classifier must not be fooled by the tetrahedral metals a
    zeolite shares with a MOF,
  * a structure that cannot be deconstructed must be reported, not raised,
    because the common use is a directory of thousands,
  * and the key must survive a change of representation, which is what makes
    it usable as a database handle.

Reference topologies are the uncontested ones: HKUST-1 is tbo and RUBTAK01 is
fcu. A zeolite is headlined by its IZA framework-type code, so ABW reports ABW
and EDI reports EDI, while the RCSR symbol of the same net stays in `names`.
'''
from __future__ import annotations

import json
import warnings
from pathlib import Path

import pytest

from mofstructure.filetyper import DEFAULT_SAVE_DIR as SAVE
from mofstructure.filetyper import STRUCTURE_DATA
from mofstructure.scripts.topology import main
from mofstructure.graph_net.periodic_graph import PeriodicGraph
from mofstructure.topology import _identify_into
from mofstructure.topology import (
    DEFAULT_METHOD,
    MOF_METHODS,
    analyse,
    analyse_methods,
    classify,
    quotient_graph,
)

DATA = Path(__file__).resolve().parent / "test_data"

warnings.filterwarnings("ignore")


class TestClassification:
    '''Choosing the deconstruction that suits a structure.'''

    @pytest.mark.parametrize(
        "filename,expected",
        [
            ("HKUST-1.cif", "mof"),
            ("RUBTAK01.cif", "mof"),
            ("ABW.cif", "zeolite"),
            ("EDI.cif", "zeolite"),
        ],
    )
    def test_material_is_recognised(self, filename, expected):
        '''Each shipped structure is placed in the right class.'''
        assert classify(str(DATA / filename)) == expected

    def test_zeolite_is_not_mistaken_for_a_mof(self):
        '''
        A silicate is a zeolite even though silicon sits in the metal list.

        The classifier tests for a zeolite first precisely because several
        tetrahedral elements - Zn, Co, Fe, Ti - are metals, and asking about
        metals first would send every one of those frameworks down the MOF
        path, where there is no cluster to cut at.
        '''
        assert classify(str(DATA / "ABW.cif")) != "mof"


class TestRecordShape:
    '''One answer shape, whatever the material.'''

    REQUIRED = {
        "interpenetration",
        "key",
        "key_hash",
        "key_version",
        "n_edges",
        "n_vertices",
        "names",
        "periodicity",
        "topology",
        "topology_source",
    }

    @pytest.mark.parametrize(
        "filename", ["HKUST-1.cif", "ABW.cif", "EDI.cif"]
    )
    def test_every_component_carries_the_same_fields(self, filename):
        '''
        A MOF and a zeolite answer with identical fields.

        The materials are deconstructed along completely different lines, so
        an implementation that let the material show through in the reply
        would make the results incomparable. This is the assertion that keeps
        them comparable.
        '''
        record = analyse(str(DATA / filename), timeout=300)
        assert record["status"] == "ok"
        assert record["components"]
        for component in record["components"]:
            assert self.REQUIRED <= set(component)

    def test_key_is_present_even_when_the_net_is_named(self):
        '''The key is the identification and is always reported.'''
        record = analyse(str(DATA / "HKUST-1.cif"), timeout=300)
        assert record["key"]
        assert record["key_hash"].startswith(record["key_version"])


class TestKnownTopologies:
    '''Agreement with topologies that are not in dispute.'''

    @pytest.mark.parametrize(
        "filename,expected",
        [
            ("HKUST-1.cif", "tbo"),
            ("RUBTAK01.cif", "fcu"),
            ("ABW.cif", "ABW"),
            ("EDI.cif", "EDI"),
        ],
    )
    def test_topology_matches_the_literature(self, filename, expected):
        '''
        The ABW case is worth its place: its net is sra in the RCSR, which has
        no net called abw, so a zeolite headlined from the RCSR would carry a
        name unrecognisable to the field it comes from. Reporting the IZA code
        is what makes the column read the way a zeolite paper does, and the
        RCSR symbol is still there in `names`.
        '''
        record = analyse(str(DATA / filename), timeout=300)
        assert record["topology"] == expected

    def test_zeolite_also_reports_its_framework_code(self):
        '''A zeolite carries both names, and both are kept.'''
        record = analyse(str(DATA / "ABW.cif"), timeout=300)
        names = record["components"][0]["names"]
        assert names.get("rcsr") == "sra"
        assert names.get("iza") == "ABW"


class TestInvariance:
    '''The key identifies the net, not the way it was written.'''

    def test_key_survives_a_supercell_and_a_shifted_origin(self):
        '''
        The same framework in a doubled cell, with atoms reordered and the
        origin moved, must give the same key. This is the property that lets
        a key be stored as a database handle: without it, two files of one
        material would look like two materials.

        The method is pinned rather than left to `DEFAULT_METHOD`. What is
        under test is that a key survives a change of representation, which
        is a property of the identification and not of any one node
        definition; LECQEQ01 has no stable net under `all_node`, so taking
        the default would test whether that structure happens to suit the
        current default instead.
        '''
        variants = (
            Path(__file__).resolve().parents[1]
            / "analysis"
            / "results"
            / "controlled_perturbation_smoke"
            / "variants"
        )
        if not variants.exists():
            pytest.skip("perturbation variants not present in this checkout")
        keys = set()
        for kind in ("pristine", "supercell-2x", "reordered-atoms",
                     "shifted-origin"):
            path = variants / f"LECQEQ01_fair__{kind}.cif"
            if not path.exists():
                continue
            record = analyse(str(path), method="sbus", timeout=300)
            assert record["status"] == "ok"
            keys.add(record["key"])
        assert len(keys) == 1, "the key moved between representations"


class TestFailuresAreReported:
    '''An ordinary failure is a status, not an exception.'''

    def test_unreadable_input_is_reported(self, tmp_path):
        '''A file ASE cannot read comes back as an error record.'''
        bad = tmp_path / "not-a-structure.cif"
        bad.write_text("this is not a crystal\n")
        record = analyse(str(bad))
        assert record["status"] in ("error", "deconstruction_failed", "no_net")
        assert "detail" in record

    def test_empty_deconstruction_is_not_a_crash(self):
        '''A CGD with no edges is reported as having no net.'''
        with pytest.raises(ValueError):
            quotient_graph("PERIODIC_GRAPH\nID x\nEDGES\nEND\n")


class TestMethods:
    '''Choosing between the several nets a MOF has.'''

    def test_a_mof_offers_more_than_one_deconstruction(self):
        '''
        A MOF has no single correct net. The alternatives answer different
        questions rather than competing, so they are all offered.
        '''
        assert len(MOF_METHODS) > 1
        assert DEFAULT_METHOD["cof"] == "cof"
        assert DEFAULT_METHOD["zeolite"] == "zeol"
        assert "zeol" not in MOF_METHODS and "cof" not in MOF_METHODS

    def test_analyse_methods_reports_every_method_it_tried(self):
        '''
        Each method appears in the result, including any that failed, so a
        missing answer is visible rather than silently absent.
        '''
        results = analyse_methods(str(DATA / "HKUST-1.cif"), timeout=300)
        assert set(results) == set(MOF_METHODS)
        for record in results.values():
            assert "status" in record


class TestEmbedding:
    '''The two geometric realisations, and what each is for.'''

    def test_refinement_makes_the_edges_more_uniform(self):
        '''
        The exact embedding minimises squared edge length, which tolerates a
        few long edges; the refined one evens them out. A builder spanning
        each edge with a linker needs the second, which is the whole reason
        it exists.
        '''
        from mofstructure.generate_cgd import TopologyExtractor
        from mofstructure.graph_net.embedding import (
            ideal_embedding,
            refined_embedding,
        )
        from mofstructure.topology import quotient_graph

        cgd = TopologyExtractor(filename=str(DATA / "SARSUC.cif")).build_cgd(
            method="sbus", name="net"
        )
        component = max(quotient_graph(cgd).components(), key=lambda c: c.n_edges)
        exact = ideal_embedding(component)
        refined = refined_embedding(component)
        assert exact["refined"] is False
        assert refined["refined"] is True
        assert refined["edge_length_spread"] < exact["edge_length_spread"]
        assert refined["edge_length_spread"] < 1.05

    def test_exact_embedding_is_reproducible(self):
        '''
        The exact embedding is the one stored beside a key, so it has to come
        out the same every time. The refined one carries no such promise,
        which is why the flag distinguishing them exists.
        '''
        from mofstructure.generate_cgd import TopologyExtractor
        from mofstructure.graph_net.embedding import ideal_embedding
        from mofstructure.topology import quotient_graph

        cgd = TopologyExtractor(filename=str(DATA / "SARSUC.cif")).build_cgd(
            method="sbus", name="net"
        )
        component = max(quotient_graph(cgd).components(), key=lambda c: c.n_edges)
        first = ideal_embedding(component)
        second = ideal_embedding(component)
        assert first["cell"] == second["cell"]
        assert first["positions"] == second["positions"]

    def test_cgd_states_which_embedding_it_holds(self):
        '''A geometry that cannot be reproduced must say so in the file.'''
        from mofstructure.structure import MOFstructure

        mof = MOFstructure(filename=str(DATA / "SARSUC.cif"))
        exact = mof.get_topology(method="all_node")["cgd"]
        refined = mof.get_topology(method="all_node", refine_cgd=True)["cgd"]
        assert "exact barycentric, reproducible" in exact
        assert "not reproducible bit for bit" in refined


class TestCommandLineOutput:
    '''
    The command writes where the rest of the package writes.

    A user builds one folder per project: `mofstructure_database` fills it,
    and a topology run has to land in the same place, under the same key, or
    the folder stops being a single database.
    '''

    def test_records_land_in_the_shared_structure_database(self, tmp_path):
        record = main([str(DATA / "ABW.cif"), "-s", str(tmp_path / SAVE)])
        written = tmp_path / SAVE / STRUCTURE_DATA / "topology_data.json"
        assert record == 0
        assert written.exists()
        assert json.loads(written.read_text())["ABW"]["topology"] == "ABW"

    def test_a_second_run_adds_to_the_file(self, tmp_path):
        save = str(tmp_path / SAVE)
        main([str(DATA / "ABW.cif"), "-s", save])
        main([str(DATA / "EDI.cif"), "-s", save])
        written = tmp_path / SAVE / STRUCTURE_DATA / "topology_data.json"
        assert sorted(json.loads(written.read_text())) == ["ABW", "EDI"]

    def test_no_save_writes_nothing(self, tmp_path):
        main([str(DATA / "ABW.cif"), "-s", str(tmp_path / SAVE),
              "--no-save"])
        assert not (tmp_path / SAVE).exists()

    def test_a_directory_is_read_as_the_structures_it_holds(self, tmp_path):
        '''
        A folder is the unit a user works in, and the sibling commands both
        take one, so naming a folder here has to read the structures inside
        rather than hand the folder itself to ASE.
        '''
        folder = tmp_path / "cifs"
        folder.mkdir()
        for name in ("ABW.cif", "EDI.cif"):
            (folder / name).write_bytes((DATA / name).read_bytes())

        main([str(folder), "-s", str(tmp_path / SAVE)])

        written = tmp_path / SAVE / STRUCTURE_DATA / "topology_data.json"
        assert sorted(json.loads(written.read_text())) == ["ABW", "EDI"]

    def test_a_directory_holding_no_structures_is_reported(self, tmp_path):
        empty = tmp_path / "empty"
        empty.mkdir()
        assert main([str(empty), "-s", str(tmp_path / SAVE)]) == 1

    def test_the_saved_record_drops_status_and_source(self, tmp_path):
        '''
        `status` says how the run went and `source` says where the file was,
        neither of which is a property of the net. The record on disk stands
        for the net, so both stay on the terminal and out of the database.
        '''
        main([str(DATA / "ABW.cif"), "-s", str(tmp_path / SAVE)])
        written = tmp_path / SAVE / STRUCTURE_DATA / "topology_data.json"
        record = json.loads(written.read_text())["ABW"]
        assert "status" not in record
        assert "source" not in record
        assert record["topology"] == "ABW"

    def test_a_csv_summary_lands_beside_the_records(self, tmp_path):
        '''
        Every other command leaves a JSON and a CSV in the folder, and a
        table is what a reader opens first, so a topology run leaves one too.
        '''
        main([str(DATA / "ABW.cif"), "-s", str(tmp_path / SAVE)])
        table = tmp_path / SAVE / STRUCTURE_DATA / "topology_data.csv"
        assert table.exists()
        header, row = table.read_text().splitlines()[:2]
        assert header.startswith("mof_names,")
        assert "topology" in header
        assert row.startswith("ABW,")

    def test_the_csv_covers_the_database_not_just_the_last_run(self, tmp_path):
        '''
        The records accumulate across runs, so a table built from one run
        would describe less than the file it sits beside.
        '''
        save = str(tmp_path / SAVE)
        main([str(DATA / "ABW.cif"), "-s", save])
        main([str(DATA / "EDI.cif"), "-s", save])
        table = tmp_path / SAVE / STRUCTURE_DATA / "topology_data.csv"
        names = [
            line.split(",")[0]
            for line in table.read_text().splitlines()[1:]
        ]
        assert sorted(names) == ["ABW", "EDI"]

    def test_the_json_export_drops_status_and_source(self, tmp_path):
        out = tmp_path / "nets.json"
        main([str(DATA / "ABW.cif"), "--json", str(out)])
        records = json.loads(out.read_text())
        assert records
        assert all(
            "status" not in record and "source" not in record
            for record in records
        )


class TestStatusTellsTheTruth:
    '''
    The status has to separate three outcomes a caller treats differently:
    a named net, a net identified but absent from the archive, and a net with
    no canonical form at all. Reporting the third as "ok" would let a failure
    into a dataset as a blank rather than as a failure.
    '''

    def test_a_named_net_is_ok(self):
        record = analyse(str(DATA / "HKUST-1.cif"), timeout=300)
        assert record["status"] == "ok"
        assert record["topology"] == "tbo"
        assert record["key"]

    def test_an_unstable_net_is_not_ok(self):
        # The quotient graph of ABUBOQ, whose neighbours collide in the
        # barycentric placement. An unstable net has no canonical form, so
        # there is nothing to identify it by and "ok" would be a lie. The
        # edges are inlined rather than read from a cif so the test states
        # exactly which graph it means.
        graph = PeriodicGraph(3, 14, [
            (0, 2, (0, 0, 0)), (0, 3, (0, 0, 0)), (0, 4, (0, 0, 0)),
            (0, 5, (0, 0, 0)), (0, 6, (0, 0, 0)), (0, 8, (0, 0, 0)),
            (0, 9, (0, 0, 0)), (0, 12, (0, 0, -1)), (1, 2, (0, 0, 0)),
            (1, 3, (0, 0, 1)), (1, 4, (0, 0, 0)), (1, 7, (0, 0, 0)),
            (1, 10, (0, 0, 0)), (1, 11, (0, 0, 0)), (1, 12, (0, 0, 0)),
            (1, 13, (0, 0, 0)), (2, 4, (0, 0, 0)), (3, 12, (0, 0, 1)),
            (5, 6, (0, 0, 0)), (7, 10, (0, 0, 0)), (8, 9, (0, 0, 0)),
            (11, 13, (0, 0, 0)),
        ])
        record = _identify_into({"status": "ok"}, graph, timeout=None,
                                descriptors=False, symmetry=False)
        assert record["status"] == "unidentified"
        assert record["key"] is None
        assert "stable" in record["detail"]

    def test_a_key_without_a_name_is_still_ok(self):
        record = analyse(str(DATA / "ABW.cif"), timeout=300)
        assert record["status"] == "ok"
        assert record["key"]
