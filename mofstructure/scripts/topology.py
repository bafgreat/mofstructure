#!/usr/bin/env python3
'''
Command-line entry point for topology identification.

Reports the net of one or more frameworks. The kind of framework is worked
out from the structure unless a deconstruction is named, so a directory
holding MOFs, COFs and zeolites together can be handed over in one go.

Records are appended to ``<save_dir>/Structure_Data/topology_data.json``,
the same file ``mofstructure_database`` fills, so one folder holds one
structure database however it was built.

A directory may be named instead of a file, in which case every structure
inside it is read, the way ``mofstructure_database`` and
``mofstructure_curate`` take a folder.

**usage:**
    mofstructure_topology HKUST-1.cif
    mofstructure_topology cif_folder
    mofstructure_topology *.cif -s MOFstructureDB
    mofstructure_topology mof.cif --method all_node
    mofstructure_topology mof.cif --all-methods
'''
from __future__ import annotations

import argparse
import json
import logging
import sys
from collections import Counter
from pathlib import Path

import mofstructure.filetyper as read_write
from mofstructure.topology import MOF_METHODS, analyse, analyse_methods

logger = logging.getLogger(__name__)

# Suffixes a named directory is expanded into. The command reads anything ASE
# reads and CGD nets besides, so the list is wider than the bare `.cif` its
# sibling commands look for, but it stays a list: a folder of structures also
# holds READMEs and JSON output, and those are not structures.
STRUCTURE_SUFFIXES = (".cif", ".cgd", ".xyz", ".pdb", ".vasp", ".res", ".gen")

# Fields dropped on the way to disk. Both describe the run rather than the
# net: `status` is what the terminal reports and what the exit code is built
# from, and `source` is the path the structure was read from, which the key
# of the record already names.
UNSAVED_FIELDS = ("source", "status")

# argparse reflows whatever it is handed, so the module docstring - written in
# reStructuredText for the API documentation - arrives on the terminal as one
# unreadable paragraph. The help text is therefore written separately.
DESCRIPTION = """\
Identify the underlying net of one or more frameworks.

The kind of framework is worked out from the structure, so a folder holding
MOFs, COFs and zeolites together can be handed over in one go. Records are
appended to <save_dir>/Structure_Data/topology_data.json, the same file
mofstructure_database fills, with a topology_data.csv summary beside it, so
one folder holds one structure database however it was built.
"""

EPILOG = """\
examples:
  mofstructure_topology HKUST-1.cif          one structure
  mofstructure_topology cif_folder           every structure in a folder
  mofstructure_topology *.cif -s MOFdb       a shell glob, saved to MOFdb
  mofstructure_topology mof.cif --method all_node
  mofstructure_topology mof.cif --all-methods

exit status:
  0  every structure was identified
  1  at least one was not; the reason is in the status column

long runs and clusters:
  Each structure runs in a worker process with a hard limit (--max-time), so
  a crash or a stall costs only that structure; it is reported as crashed or
  timeout. Results are checkpointed under
  <save_dir>/Structure_Data/_progress/ as they finish, and repeating the
  command continues a killed run. Use -j for several workers and
  --shard slurm/N on a cluster array, then mofstructure_merge <save_dir>.
"""


def _database_path(save_dir: str) -> Path:
    '''
    Build the path records are appended to, laid out the way every other
    mofstructure command lays out its output.

    **parameters:**
        - save_dir: str
            Root of the structure database.

    **returns:**
        pathlib.Path
            `<save_dir>/Structure_Data/topology_data.json`, its parent
            created when missing.
    '''
    folder = Path(save_dir) / read_write.STRUCTURE_DATA
    folder.mkdir(parents=True, exist_ok=True)
    return folder / "topology_data.json"


def _structure_files(folder: Path) -> list[Path]:
    '''
    List the structures a named directory holds.

    The search does not descend, matching the way `mofstructure_database`
    and `mofstructure_curate` read a folder. Hidden files are skipped, which
    on macOS also drops the `._` companions an archive leaves behind and
    that ASE cannot read.

    **parameters:**
        - folder: pathlib.Path
            Directory to list.

    **returns:**
        list of pathlib.Path
            Structure files, in name order.
    '''
    return sorted(
        item
        for item in folder.iterdir()
        if item.is_file()
        and not item.name.startswith(".")
        and item.suffix.lower() in STRUCTURE_SUFFIXES
    )


def _for_disk(record: dict) -> dict:
    '''
    Drop the fields that describe the run rather than the net.

    **parameters:**
        - record: python dictionary
            Result of `analyse`, as printed.

    **returns:**
        python dictionary
            A copy without `UNSAVED_FIELDS`.
    '''
    return {
        key: value for key, value in record.items()
        if key not in UNSAVED_FIELDS
    }


def _csv_row(record: dict) -> dict:
    '''
    Flatten one record into the columns a summary table wants.

    A record nests one entry per component, which a table cannot hold, so
    the largest component is reported: a framework in several pieces is the
    framework proper plus interpenetrating copies of it, and the largest is
    the framework. This is the choice `MOFstructure.get_topology` makes, so
    both writers of the folder agree on what the row means.

    **parameters:**
        - record: python dictionary
            Result of `analyse`.

    **returns:**
        python dictionary
            One flat row, without the fields `_for_disk` drops.
    '''
    row = {
        "material": record.get("material"),
        "method": record.get("method"),
        "topology": record.get("topology"),
        "topology_source": None,
        "dimension": None,
        "key_hash": record.get("key_hash"),
        "n_components": record.get("n_components"),
        "detail": record.get("detail"),
    }
    components = record.get("components") or []
    if components:
        largest = max(components, key=lambda c: c.get("n_edges") or 0)
        row["topology_source"] = largest.get("topology_source")
        row["dimension"] = largest.get("periodicity")
        row["key_hash"] = row["key_hash"] or largest.get("key_hash")
    return row


def _write_csv(records: dict, destination: Path) -> Path:
    '''
    Write the summary table beside the records it summarises.

    **parameters:**
        - records: mapping
            Every record in the database, keyed by structure name.

        - destination: pathlib.Path
            The JSON the records were written to.

    **returns:**
        pathlib.Path
            The CSV written, named after the JSON.
    '''
    table = destination.with_suffix(".csv")
    rows = {
        name: _csv_row(record)
        for name, record in records.items()
        if isinstance(record, dict) and record
    }
    read_write.summary_frame(rows).to_csv(table)
    return table


def _row(record: dict, label: str) -> str:
    '''
    Render one record as a single line.

    **parameters:**
        - record: python dictionary
            Result of `analyse`.

        - label: str
            Name to show for the structure.

    **returns:**
        str
    '''
    if record["status"] != "ok":
        return (
            f"{label:<28} {record['status']:<22} "
            f"{str(record.get('detail', ''))[:60]}"
        )
    parts = []
    for component in record["components"]:
        name = component.get("topology") or "unnamed"
        source = component.get("topology_source")
        parts.append(f"{name}({source})" if source else name)
    return (
        f"{label:<28} {'ok':<22} "
        f"{record.get('material') or '?':<8} {record['method']:<14} "
        f"{' + '.join(parts)}"
    )


def main(argv: list[str] = None) -> int:
    '''
    Identify the topology of every structure named on the command line.

    **parameters:**
        - argv: list of str, optional
            Arguments; `sys.argv[1:]` when omitted.

    **returns:**
        int
            0 when every structure was identified, 1 otherwise.
    '''
    parser = argparse.ArgumentParser(
        description=DESCRIPTION,
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "structures",
        nargs="+",
        type=Path,
        help=(
            "structure files, or a directory, in which case every "
            f"{'/'.join(STRUCTURE_SUFFIXES)} file it holds is read"
        ),
    )
    parser.add_argument(
        "--method",
        default="auto",
        help=(
            "deconstruction to use; 'auto' works it out from the structure. "
            f"A MOF also accepts {', '.join(MOF_METHODS)}"
        ),
    )
    parser.add_argument(
        "--all-methods",
        action="store_true",
        help=(
            "report every alternative net a MOF has; a COF or a "
            "zeolite has one"
        ),
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=300,
        help="seconds allowed per structure, 0 for no limit",
    )
    parser.add_argument(
        "--descriptors",
        action="store_true",
        help="also record the coordination sequence and density of the net",
    )
    parser.add_argument(
        "--symmetry",
        action="store_true",
        help="also record the maximal space group and intrinsic chirality",
    )
    parser.add_argument(
        "-s",
        "--save_dir",
        default=read_write.DEFAULT_SAVE_DIR,
        help=(
            "directory to save output files, laid out as the other "
            "mofstructure commands lay it out: records are appended to "
            "<save_dir>/Structure_Data/topology_data.json, keyed by "
            "structure "
            "name, so a topology run and a database run fill the same file"
        ),
    )
    parser.add_argument(
        "--json",
        type=Path,
        help="write the records to this file instead of the save directory",
    )
    parser.add_argument(
        "--no-save",
        action="store_true",
        help="print only, writing nothing",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="report each deconstruction as it is built",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="print the closing summary only, without a line per structure",
    )
    from mofstructure.batch_tasks import add_batch_arguments

    add_batch_arguments(parser, max_time=1800)
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s - %(message)s",
        force=True,
    )

    options = {
        "timeout": args.timeout or None,
        "descriptors": args.descriptors,
        "symmetry": args.symmetry,
    }

    records = []
    collected = {}
    failures = 0

    # Directories are resolved before anything is read, so a folder holding
    # nothing readable is reported as that rather than handed to ASE, which
    # would call it a corrupt trajectory.
    targets: list[Path] = []
    for path in args.structures:
        if path.is_dir():
            found = _structure_files(path)
            if not found:
                print(f"{str(path):<28} no structures found")
                failures += 1
            targets.extend(found)
        elif path.exists():
            targets.append(path)
        else:
            print(f"{str(path):<28} missing")
            failures += 1

    # Every structure runs in a worker process with a hard time limit, so a
    # crash or a stall in one structure cannot end the run, and results are
    # checkpointed as they finish so a killed job resumes where it stopped.
    import tempfile
    from types import SimpleNamespace

    from mofstructure.batch import (parse_shard, run_batch, select_shard,
                                    shard_suffix)
    from mofstructure.batch_tasks import consolidate, done_names

    label = "all-methods" if args.all_methods else args.method
    command = f"topology.{label}"
    if args.no_save or args.json:
        scratch = tempfile.mkdtemp(prefix="mofstructure_topology_")
        structure_db = Path(scratch) / read_write.STRUCTURE_DATA
    else:
        scratch = None
        structure_db = _database_path(args.save_dir).parent
    shard = parse_shard(args.shard)
    items = select_shard([(path.stem, str(path)) for path in targets], shard)
    skip = done_names(structure_db, command, retry_failed=args.retry_failed)
    total = len(items)
    counter = len(str(total)) * 2 + 2 if total > 1 else 0
    if skip and not args.quiet:
        print(f"{len(skip & {n for n, _ in items})} structures already "
              f"recorded by an earlier run are skipped")
    if items and not args.quiet:
        print(
            f"{'':<{counter}}{'structure':<28} {'status':<22} "
            f"{'kind':<8} {'method':<14} net", flush=True
        )

    def unpack(batch_record):
        '''The analyse record(s) a checkpoint line stands for.'''
        name = batch_record["name"]
        if batch_record["status"] != "ok":
            failed = {"status": batch_record["status"],
                      "detail": batch_record.get("detail"),
                      "method": label, "components": []}
            return [(name, failed)]
        data = batch_record["data"]
        if args.all_methods:
            return [(f"{name}:{method}", sub) for method, sub in data.items()]
        return [(name, data)]

    def show(batch_record, done, _total):
        nonlocal failures
        stamp = f"{done:>{len(str(total))}}/{total} " if counter else ""
        for key, record in unpack(batch_record):
            records.append(record)
            collected[key] = record
            failures += record["status"] != "ok"
            if not args.quiet:
                shown = f"{Path(batch_record['file']).name}" + (
                    f":{key.split(':', 1)[1]}" if ":" in key else "")
                print(stamp + _row(record, shown), flush=True)

    run_batch(
        items,
        "mofstructure.batch_tasks:topology_task",
        structure_db / "_progress" / f"{command}{shard_suffix(shard)}.jsonl",
        kwargs={"method": args.method, "all_methods": args.all_methods,
                "timeout": options["timeout"],
                "descriptors": options["descriptors"],
                "symmetry": options["symmetry"]},
        workers=args.workers,
        timeout=args.max_time or None,
        memory_gb=args.memory_limit,
        recycle=args.recycle,
        skip=skip,
        on_record=show,
    )

    if not args.no_save and records:
        if args.json:
            args.json.parent.mkdir(parents=True, exist_ok=True)
            args.json.write_text(
                json.dumps(
                    [_for_disk(record) for record in records],
                    indent=1,
                    default=str,
                )
            )
            destination = args.json
            table = _write_csv(
                {
                    name: _for_disk(record)
                    for name, record in collected.items()
                },
                destination,
            )
        elif shard:
            destination = structure_db / "_progress"
            table = None
            print(f"\nshard written to {destination}; join all shards with "
                  f"mofstructure_merge {args.save_dir}")
        else:
            consolidate(args.save_dir)
            destination = _database_path(args.save_dir)
            # Built from the merged file rather than from this run, so the
            # table covers the database and not just the last command.
            table = _write_csv(
                read_write.load_data(str(destination)), destination
            )
        print(f"\nwrote {len(records)} records to {destination}")
        if table:
            print(f"summary table {table}")
    if scratch:
        import shutil

        shutil.rmtree(scratch, ignore_errors=True)

    if records:
        tally = Counter(record["status"] for record in records)
        print(
            f"\n{len(records)} identified: "
            + ", ".join(
                f"{count} {status}" for status, count in tally.most_common()
            )
        )

    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
