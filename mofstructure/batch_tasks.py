#!/usr/bin/env python3
'''
Per-structure tasks run by the batch commands, and the step that joins their
checkpoint files into the usual structure database.

Each task takes the path of one structure and returns json-serialisable data.
They are run by `mofstructure.batch.run_batch` in worker processes, so an
exception is recorded for that structure and a crash or a stall costs only
that structure.

**database layout:**
    The commands write json-lines checkpoints to

        <save_dir>/Structure_Data/_progress/<command>[.part-iiii-of-nnnn].jsonl

    with <command> one of "database", "topology", "porosity" and "oms", and
    `consolidate` turns them into the files the commands have always written:

        Structure_Data/sbus_and_linkers.json         database
        Structure_Data/ligands_data.json             database
        Structure_Data/porosity_data.json            database, porosity
        Structure_Data/structure_oms_and_general_info.json   database --oms, oms
        Structure_Data/topology_data.json            database, topology
        Structure_Data/fingerprint_data.json         database
        Structure_Data/run_status.csv                every command

    each json with a csv summary beside it. run_status.csv has one row per
    structure and command with the status (ok, error, timeout, crashed), the
    time taken and the reason for a failure, so a structure that did not
    finish is visible instead of missing.
'''
from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd

import mofstructure.filetyper as read_write
from mofstructure.batch import CHECKPOINT_DIR, read_checkpoints

DATABASE_FILES = {
    "sbu": "sbus_and_linkers.json",
    "ligand": "ligands_data.json",
    "porosity": "porosity_data.json",
    "oms": "structure_oms_and_general_info.json",
    "topology": "topology_data.json",
    "fingerprint": "fingerprint_data.json",
}

# Fields of a topology record that describe the run rather than the net.
TOPOLOGY_UNSAVED = ("source", "status")


def _plain(data):
    '''Round-trip through json so numpy scalars and similar become plain.'''
    return json.loads(json.dumps(read_write.convert_numpy_types(data),
                                 default=str))


def database_name(path: str) -> str:
    '''Structure name used by the database commands: text before first dot.'''
    return os.path.basename(path).split(".")[0]


# ---------------------------------------------------------------------------
# Tasks
# ---------------------------------------------------------------------------

def topology_task(path, method="auto", all_methods=False, timeout=300,
                  descriptors=False, symmetry=False):
    '''Topology record of one structure, as `mofstructure_topology` reports.'''
    from mofstructure.topology import analyse, analyse_methods

    options = {"timeout": timeout or None, "descriptors": descriptors,
               "symmetry": symmetry}
    if all_methods:
        return _plain(analyse_methods(str(path), **options))
    return _plain(analyse(str(path), method=method, **options))


def porosity_task(path, probe_radius=1.86, number_of_steps=10000,
                  rad_file=None, high_accuracy=True, timeout=1800):
    '''Pore geometry of one guest-free structure.'''
    from mofstructure import structure

    mof = structure.MOFstructure(filename=str(path))
    return _plain(mof.get_porosity(
        probe_radius=probe_radius, number_of_steps=number_of_steps,
        rad_file=rad_file, high_accuracy=high_accuracy, timeout=timeout))


def oms_task(path, max_atoms=5000):
    '''Open metal sites and general information of one structure.'''
    from mofstructure import structure

    mof = structure.MOFstructure(filename=str(path))
    if max_atoms and len(mof.ase_atoms) > max_atoms:
        raise ValueError(f"skipped: {len(mof.ase_atoms)} atoms is more than "
                         f"--max-atoms {max_atoms}")
    return _plain(mof.get_oms())


def database_task(path, xyz_path, oms=False, topology=True,
                  topology_method="auto", topology_timeout=300,
                  porosity_timeout=1800, probe_radius=1.86,
                  high_accuracy=True, oms_max_atoms=5000, only=None):
    '''
    Every analysis `mofstructure_database` records for one structure.

    Each analysis is attempted on its own, so a failure in one leaves the
    others in place; the failures are listed under "errors". `only` limits
    the run to the named analyses, which is how a database written by an
    earlier version is filled in without repeating what it already holds.
    '''
    from mofstructure.scripts import collect_data

    name = database_name(path)
    mof = collect_data.structure.MOFstructure(filename=str(path))
    out, errors = {}, {}
    wanted = set(only) if only else None

    def attempt(key, function):
        if wanted is not None and key not in wanted:
            return
        try:
            out[key] = _plain(function())
        except Exception as exc:  # noqa: BLE001 - recorded per analysis
            out[key] = None
            errors[key] = f"{type(exc).__name__}: {exc}"[:300]

    if topology:
        def net():
            record = mof.get_topology(method=topology_method,
                                      timeout=topology_timeout)
            return {k: v for k, v in record.items() if k != "status"}
        attempt("topology", net)

    attempt("fingerprint", mof.get_ligand_cluster_fingerprint)

    def sbu():
        found = mof.get_sbu()
        return None if found is None else collect_data.collect_sbus(
            *found, name, xyz_path)
    attempt("sbu", sbu)

    def ligand():
        found = mof.get_ligands()
        return None if found is None else collect_data.collect_ligand(
            found[1], name, xyz_path)
    attempt("ligand", ligand)

    attempt("porosity", lambda: mof.get_porosity(
        probe_radius=probe_radius, timeout=porosity_timeout,
        high_accuracy=high_accuracy))

    if oms and (wanted is None or "oms" in wanted):
        if len(mof.ase_atoms) > oms_max_atoms:
            out["oms"] = None
            errors["oms"] = f"skipped: more than {oms_max_atoms} atoms"
        else:
            attempt("oms", mof.get_oms)

    if errors:
        out["errors"] = errors
    return out


# ---------------------------------------------------------------------------
# Joining checkpoints into the database
# ---------------------------------------------------------------------------

def checkpoint_files(structure_db: Path, command: str) -> list[Path]:
    '''
    Every checkpoint part a command has written, whole run or shards.

    A command name ending in ".*" matches every variant, which is how the
    topology checkpoints, written one file per deconstruction method, are
    collected together.
    '''
    folder = Path(structure_db) / CHECKPOINT_DIR
    if command.endswith(".*"):
        stem = command[:-2]
        return sorted(folder.glob(f"{stem}.jsonl")) + \
            sorted(folder.glob(f"{stem}.*.jsonl"))
    return sorted(folder.glob(f"{command}.jsonl")) + \
        sorted(folder.glob(f"{command}.part-*.jsonl"))


def done_names(structure_db: Path, command: str,
               retry_failed: bool = False) -> set[str]:
    '''Names a command has already recorded, for resuming a run.'''
    records = read_checkpoints(checkpoint_files(structure_db, command))
    if retry_failed:
        return {n for n, r in records.items() if r.get("status") == "ok"}
    return set(records)


def _load(path: Path) -> dict:
    if path.exists() and path.stat().st_size:
        return read_write.load_data(str(path))
    return {}


def _write(records: dict, path: Path, drop=()) -> None:
    read_write.write_json_atomic(records, str(path))
    frame = read_write.summary_frame(records, drop=drop)
    if not frame.empty:
        frame.to_csv(path.with_suffix(".csv"))


def consolidate(save_dir: str) -> dict[str, int]:
    '''
    Join every checkpoint under `save_dir` into the structure database files.

    Records already in the json files, for example from a run made with an
    older version, are kept, and a checkpoint record replaces the stored one
    for the same structure. Safe to run at any time, including while a run is
    in progress: it only reads the checkpoints.

    **returns:**
        python dictionary
            Number of structures written per database file.
    '''
    from mofstructure.porosity import empty_porosity_record

    structure_db = Path(save_dir) / read_write.STRUCTURE_DATA
    structure_db.mkdir(parents=True, exist_ok=True)
    status_rows = []
    tables: dict[str, dict] = {}

    def table(key):
        if key not in tables:
            tables[key] = _load(structure_db / DATABASE_FILES[key])
        return tables[key]

    def note(command, record):
        status_rows.append({
            "name": record["name"], "command": command,
            "status": record.get("status"),
            "elapsed_s": record.get("elapsed_s"),
            "detail": record.get("detail"),
            "analysis_errors": json.dumps(
                (record.get("data") or {}).get("errors"))
            if command == "database" and isinstance(record.get("data"), dict)
            else None,
            "file": record.get("file"),
        })

    for record in read_checkpoints(checkpoint_files(structure_db,
                                                    "database")).values():
        note("database", record)
        name = record["name"]
        if record.get("status") == "ok":
            for key, value in record["data"].items():
                if key in DATABASE_FILES:
                    table(key)[name] = value
        else:
            table("porosity")[name] = empty_porosity_record(record["status"])

    topology_parts = checkpoint_files(structure_db, "topology.*")
    for record in read_checkpoints(topology_parts).values():
        note("topology", record)
        if record.get("status") != "ok":
            continue
        data = record["data"]
        if data and all(isinstance(v, dict) and "method" in v
                        for v in data.values()):
            for method, sub in data.items():
                table("topology")[f"{record['name']}:{method}"] = {
                    k: v for k, v in sub.items() if k not in TOPOLOGY_UNSAVED}
        else:
            table("topology")[record["name"]] = {
                k: v for k, v in data.items() if k not in TOPOLOGY_UNSAVED}

    for record in read_checkpoints(checkpoint_files(structure_db,
                                                    "porosity")).values():
        note("porosity", record)
        if record.get("status") == "ok":
            table("porosity")[record["name"]] = record["data"]
        else:
            table("porosity")[record["name"]] = empty_porosity_record(
                record["status"])

    for record in read_checkpoints(checkpoint_files(structure_db,
                                                    "oms")).values():
        note("oms", record)
        if record.get("status") == "ok":
            table("oms")[record["name"]] = record["data"]

    written = {}
    for key, records in tables.items():
        path = structure_db / DATABASE_FILES[key]
        drop = ["cgd"] if key == "topology" else ()
        if key == "fingerprint":
            read_write.write_json_atomic(records, str(path))
            frame = read_write.summary_frame(records)
            if not frame.empty:
                columns = [c for c in ("fingerprint_hash", "cluster_units")
                           if c in frame]
                frame[columns].to_csv(path.with_suffix(".csv"))
        else:
            _write(records, path, drop=drop)
        written[DATABASE_FILES[key]] = len(records)

    if status_rows:
        pd.DataFrame(status_rows).to_csv(structure_db / "run_status.csv",
                                         index=False)
        written["run_status.csv"] = len(status_rows)
    return written


# ---------------------------------------------------------------------------
# Command-line helpers shared by the batch commands
# ---------------------------------------------------------------------------

BATCH_EPILOG = """\
long runs and clusters:
  Every structure runs in a separate worker process with a hard time limit,
  so a crash or a stall costs only that structure. Results are appended to
  <save_dir>/Structure_Data/_progress/*.jsonl as each structure finishes, and
  repeating the same command continues where a killed job stopped. The
  database files and run_status.csv are rebuilt at the end of every run and
  can be rebuilt at any time with mofstructure_merge <save_dir>.

  Split a folder across a cluster array with --shard slurm/N (or i/N); every
  task writes its own part, and mofstructure_merge joins them.
"""


def add_batch_arguments(parser, max_time=7200):
    '''Add the options every batch command shares.'''
    group = parser.add_argument_group("batch execution")
    group.add_argument("-j", "--workers", type=int, default=1,
                       help="worker processes (default 1)")
    group.add_argument("--max-time", type=float, default=max_time,
                       help="hard wall-clock limit per structure in seconds; "
                            "the worker is killed and the structure recorded "
                            f"as a timeout (default {max_time:g}, 0 for none)")
    group.add_argument("--memory-limit", type=float, default=None,
                       help="memory limit per worker in GB (Linux)")
    group.add_argument("--shard", type=str, default=None,
                       help="process one part of the folder, given as i/N "
                            "(0-based) or slurm/N to read "
                            "SLURM_ARRAY_TASK_ID")
    group.add_argument("--retry-failed", action="store_true",
                       help="run again structures recorded as error, "
                            "timeout or crashed")
    group.add_argument("--recycle", type=int, default=50,
                       help="structures per worker before it is replaced "
                            "(default 50)")
    group.add_argument("--merge-every", type=float, default=60,
                       help="rebuild the json and csv files from the "
                            "checkpoints every this many minutes during the "
                            "run (default 60, 0 only at the end); ignored "
                            "with --shard")
    return parser


def periodic_merge(save_dir, minutes, on_record=None):
    '''
    Wrap a progress callback so that the database files are rebuilt from the
    checkpoints at most every `minutes` minutes while a run is going.

    The checkpoints grow by one line per structure throughout the run; the
    json files are single objects that must be rewritten whole, so they are
    refreshed on a timer rather than after every structure.
    '''
    import time

    last = [time.monotonic()]

    def callback(record, done, total):
        if on_record:
            on_record(record, done, total)
        if minutes and done < total and \
                time.monotonic() - last[0] >= minutes * 60:
            written = consolidate(save_dir)
            last[0] = time.monotonic()
            print(f"database files refreshed after {done}/{total} "
                  f"structures ({sum(written.values())} records)", flush=True)

    return callback


def run_command(command, items, task, kwargs, save_dir, args,
                on_record=None, extra_skip=()):
    '''
    Run a batch command with the shared options and rebuild the database.

    **returns:**
        python dictionary
            Count of records by status for this run.
    '''
    from mofstructure.batch import (parse_shard, progress_printer, run_batch,
                                    select_shard, shard_suffix)

    structure_db = Path(save_dir) / read_write.STRUCTURE_DATA
    shard = parse_shard(args.shard)
    selected = select_shard(items, shard)
    skip = done_names(structure_db, command, retry_failed=args.retry_failed)
    skip |= set(extra_skip)
    checkpoint = structure_db / CHECKPOINT_DIR / \
        f"{command}{shard_suffix(shard)}.jsonl"
    todo = sum(1 for name, _ in selected if name not in skip)
    print(f"{command}: {len(selected)} structures"
          f"{' in this shard' if shard else ''}, "
          f"{len(selected) - todo} already recorded, {todo} to run, "
          f"{args.workers} worker(s)", flush=True)
    callback = on_record or progress_printer()
    minutes = getattr(args, "merge_every", 0)
    if minutes and not shard:
        callback = periodic_merge(save_dir, minutes, callback)
    tally = run_batch(
        selected, task, checkpoint, kwargs=kwargs, workers=args.workers,
        timeout=args.max_time or None, memory_gb=args.memory_limit,
        recycle=args.recycle, skip=skip, on_record=callback,
    )
    print(f"{command}: this run {tally or 'nothing to do'}", flush=True)
    if shard:
        # Array tasks finish at nearly the same time, and each rebuilding
        # the same large files would only repeat the work, so the parts are
        # joined once, after the last task, by mofstructure_merge.
        print(f"shard written to {checkpoint}; join all shards with "
              f"mofstructure_merge {save_dir}", flush=True)
        return tally
    written = consolidate(save_dir)
    for name, count in written.items():
        print(f"  {structure_db / name}: {count} records", flush=True)
    return tally
