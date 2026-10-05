"""
The batch runner must finish a run whatever individual structures do.

Each misbehaving task stands for a failure seen in long database runs: an
exception, an abort or a segmentation fault in compiled code, and a structure
that never finishes. None of them may stop the run, every structure must be
recorded with what happened to it, and a repeated run must continue rather
than start again.
"""
import json
import os

import pytest

from mofstructure import filetyper
from mofstructure.batch import (parse_shard, read_checkpoints, run_batch,
                                select_shard)

TASK = "tests._batch_helpers:behave"


def _items(*modes):
    return [(mode, f"/nowhere/{mode}") for mode in modes]


def test_failures_are_recorded_and_the_run_continues(tmp_path):
    checkpoint = tmp_path / "run.jsonl"
    items = _items("ok-1", "raise-1", "abort-1", "segfault-1", "hang-1",
                   "ok-2", "ok-3")
    tally = run_batch(items, TASK, checkpoint, workers=2, timeout=5,
                      kwargs={"flag": 1})
    records = read_checkpoints([checkpoint])
    assert set(records) == {name for name, _ in items}
    assert records["ok-1"]["status"] == "ok"
    assert records["ok-1"]["data"] == {"mode": "ok-1", "flag": 1}
    assert records["ok-3"]["status"] == "ok"
    assert records["raise-1"]["status"] == "error"
    assert "deliberate failure" in records["raise-1"]["detail"]
    assert records["abort-1"]["status"] == "crashed"
    assert "SIGABRT" in records["abort-1"]["detail"]
    assert records["segfault-1"]["status"] == "crashed"
    assert records["hang-1"]["status"] == "timeout"
    assert sum(tally.values()) == len(items)


def test_a_repeated_run_continues_instead_of_starting_again(tmp_path):
    checkpoint = tmp_path / "run.jsonl"
    run_batch(_items("ok-1", "raise-1"), TASK, checkpoint, timeout=5)
    done = set(read_checkpoints([checkpoint]))
    tally = run_batch(_items("ok-1", "raise-1", "ok-2"), TASK, checkpoint,
                      timeout=5, skip=done)
    assert tally == {"ok": 1}
    assert len(checkpoint.read_text().splitlines()) == 3


def test_a_line_cut_short_by_a_killed_job_is_ignored(tmp_path):
    checkpoint = tmp_path / "run.jsonl"
    run_batch(_items("ok-1"), TASK, checkpoint, timeout=5)
    with open(checkpoint, "a", encoding="utf-8") as handle:
        handle.write('{"name": "ok-2", "status": "o')
    assert set(read_checkpoints([checkpoint])) == {"ok-1"}


def test_shards_cover_every_structure_once():
    items = _items(*(f"ok-{i}" for i in range(23)))
    parts = [select_shard(items, (i, 4)) for i in range(4)]
    names = [name for part in parts for name, _ in part]
    assert sorted(names) == sorted(name for name, _ in items)
    assert len(names) == len(set(names))


def test_shard_index_can_come_from_the_scheduler(monkeypatch):
    monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "3")
    assert parse_shard("slurm/8") == (3, 8)
    with pytest.raises(ValueError):
        parse_shard("8/8")


def test_an_interrupted_json_write_leaves_the_old_file(tmp_path):
    target = tmp_path / "data.json"
    filetyper.append_json({"a": 1}, str(target))

    class Broken(json.JSONEncoder):
        def default(self, o):
            raise RuntimeError("killed mid write")

    with pytest.raises(RuntimeError):
        filetyper.write_json_atomic({"a": 2, "b": object()}, str(target),
                                    encoder=Broken)
    assert json.loads(target.read_text()) == {"a": 1}
    assert not [p for p in os.listdir(tmp_path) if p.startswith(".tmp_")]


def test_a_corrupt_database_file_is_reported_not_overwritten(tmp_path):
    target = tmp_path / "data.json"
    target.write_text('{"a": 1, "b":')
    with pytest.raises(ValueError, match="not valid json"):
        filetyper.append_json({"c": 3}, str(target))
    assert target.read_text() == '{"a": 1, "b":'


def test_a_timeout_kills_the_children_of_the_task(tmp_path):
    '''
    Porosity runs zeo++ in a child of the worker. Killing only the worker on
    a timeout left that child running with nothing to stop it, so a long run
    filled up with orphaned zeo++ processes until the machine stalled.
    '''
    import time

    pid_file = tmp_path / "child.pid"
    tally = run_batch(_items("child-1"), TASK, tmp_path / "run.jsonl",
                      timeout=5, kwargs={"pid_file": str(pid_file)})
    assert tally == {"timeout": 1}
    pid = int(pid_file.read_text())
    for _ in range(50):
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return
        time.sleep(0.1)
    os.kill(pid, 9)
    pytest.fail("the task's child outlived its worker")
