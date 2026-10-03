#!/usr/bin/env python3
'''
Command-line entry point `mofstructure_merge`.

Joins the checkpoint files written by `mofstructure_database`,
`mofstructure_topology` and `mofstructure_porosity`, including the parts of a
run split across a cluster array with --shard, into the json and csv files
under <save_dir>/Structure_Data, and writes run_status.csv with one row per
structure and command.

**usage:**
    mofstructure_merge MOFstructureDB
'''
import argparse
import sys

from mofstructure.batch_tasks import consolidate


def main(argv=None):
    '''Rebuild the database files from the checkpoints in a save directory.'''
    parser = argparse.ArgumentParser(
        description="Join batch checkpoints, including the shards of an "
                    "array job, into the structure database files.")
    parser.add_argument("save_dir", help="database directory to rebuild")
    args = parser.parse_args(argv)
    written = consolidate(args.save_dir)
    if not written:
        print(f"no checkpoints found under {args.save_dir}")
        return 1
    for name, count in written.items():
        print(f"{name:<45} {count} records")
    return 0


if __name__ == "__main__":
    sys.exit(main())
