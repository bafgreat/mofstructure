#!/usr/bin/python
'''
Command line entry point for `mofstructure_oms`.

Reports the open metal sites of every structure in a folder, along with the
general information about each framework, and collects them into
`<save_dir>/Structure_Data/structure_oms_and_general_info.json` and its csv,
the file `mofstructure_database --oms` fills. Every structure runs in a worker
process with a hard time limit, results are checkpointed as they finish and
repeating the command continues an interrupted run (see
`mofstructure.batch`).
'''
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
import os
import sys
import argparse
import mofstructure.filetyper as read_write


def compile_data(cif_files, result_folder, verbose=False, max_atoms=5000,
                 workers=1, max_time=3600, memory_limit=None, shard=None,
                 retry_failed=False, recycle=50, merge_every=60):
    '''
    Remove guests and record the open metal sites of a list of structures.

    **parameters:**
        - cif_files: list of str

        - result_folder: str
            Directory the database is written to.

        - verbose: bool

        - max_atoms: int
            Structures with more atoms are skipped, since the analysis can
            exhaust memory on very large cells; the skip is recorded.

        - workers, max_time, memory_limit, shard, retry_failed, recycle,
          merge_every:
            batch settings, see `mofstructure.batch_tasks.add_batch_arguments`.
    '''
    from types import SimpleNamespace
    from mofstructure.batch_tasks import database_name, run_command

    os.makedirs(os.path.join(result_folder, read_write.STRUCTURE_DATA),
                exist_ok=True)
    items = [(database_name(f), f) for f in cif_files]
    options = SimpleNamespace(workers=workers, max_time=max_time,
                              memory_limit=memory_limit, shard=shard,
                              retry_failed=retry_failed, recycle=recycle,
                              merge_every=merge_every)
    run_command('oms', items, 'mofstructure.batch_tasks:oms_task',
                {'max_atoms': max_atoms}, result_folder, options)
    if verbose:
        print(f"Saved results to {result_folder}")


def main():
    '''
    mofstructure command line interface for computing open metal sites
    '''
    from mofstructure.batch_tasks import BATCH_EPILOG, add_batch_arguments

    parser = argparse.ArgumentParser(
        description='Record the open metal sites of a folder of structures.',
        epilog=BATCH_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('cif_folder', type=str,
                        help='folder of cif files, or one cif file')
    parser.add_argument('-s', '--save_dir', type=str,
                        default=read_write.DEFAULT_SAVE_DIR,
                        help='directory to save output files')
    parser.add_argument('--max-atoms', type=int, default=5000,
                        help='skip structures with more atoms '
                             '(default 5000, 0 for no limit)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='print verbose output')
    add_batch_arguments(parser, max_time=3600)
    args = parser.parse_args()
    if os.path.isdir(args.cif_folder):
        cif_files = [os.path.join(args.cif_folder, f)
                     for f in os.listdir(args.cif_folder)
                     if f.endswith('.cif') and not f.startswith('.')]
    else:
        cif_files = [args.cif_folder]
    compile_data(cif_files, args.save_dir, args.verbose, args.max_atoms,
                 workers=args.workers, max_time=args.max_time,
                 memory_limit=args.memory_limit, shard=args.shard,
                 retry_failed=args.retry_failed, recycle=args.recycle,
                 merge_every=args.merge_every)


if __name__ == '__main__':
    sys.exit(main())
