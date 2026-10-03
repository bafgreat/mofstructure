#!/usr/bin/python
'''
Command line entry point for `mofstructure_porosity`.

Computes the pore geometry of every structure in a folder and collects the
results into one table. Guests are removed before measuring. zeo++ runs under
a timeout, and a structure it cannot handle is recorded with its reason rather
than dropped, so the table has one row per input whatever happens.
'''
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
import os
import argparse
import mofstructure.filetyper as read_write
from mofstructure.porosity import DEFAULT_TIMEOUT as porosity_timeout


def compile_data(cif_files,
                 result_folder,
                 probe_radius=1.86,
                 number_of_steps=10000,
                 rad_file=None,
                 verbose=False,
                 timeout=porosity_timeout,
                 high_accuracy=True,
                 workers=1,
                 max_time=None,
                 memory_limit=None,
                 shard=None,
                 retry_failed=False,
                 recycle=50, merge_every=60):
    '''
    Remove guests and compute the pore geometry of every structure in a list.

    Every structure runs in a separate worker process (see
    `mofstructure.batch`), so a crash or a stall in zeo++ or in the guest
    removal costs only that structure. Results are checkpointed as they
    finish and `Structure_Data/porosity_data.json` and its csv are rebuilt at
    the end; a structure that failed has a row of missing values and its
    reason in `porosity_status`. Repeating the call continues an interrupted
    run.

    **parameters:**
        - cif_files: list of str

        - result_folder: str
            Directory the database is written to.

        - probe_radius, number_of_steps, rad_file, high_accuracy:
            zeo++ settings, see `MOFstructure.get_porosity`.

        - timeout: float
            Seconds allowed for zeo++ per structure.

        - workers, max_time, memory_limit, shard, retry_failed, recycle:
            batch settings, see `mofstructure.batch_tasks.add_batch_arguments`.
            max_time defaults to the zeo++ timeout plus ten minutes for the
            guest removal.
    '''
    from types import SimpleNamespace
    from mofstructure.batch_tasks import database_name, run_command

    os.makedirs(os.path.join(result_folder, read_write.STRUCTURE_DATA),
                exist_ok=True)
    items = [(database_name(f), f) for f in cif_files]
    options = SimpleNamespace(
        workers=workers,
        max_time=max_time if max_time is not None else timeout + 600,
        memory_limit=memory_limit, shard=shard,
        retry_failed=retry_failed, recycle=recycle,
                              merge_every=merge_every)
    kwargs = dict(probe_radius=probe_radius, number_of_steps=number_of_steps,
                  rad_file=rad_file, high_accuracy=high_accuracy,
                  timeout=timeout)
    run_command('porosity', items, 'mofstructure.batch_tasks:porosity_task',
                kwargs, result_folder, options)
    if verbose:
        print(f"Saved results to {result_folder}")


def main():
    '''
    mofstructure command line interface to compute the porosity of any
    periodic system.
    '''
    from mofstructure.batch_tasks import BATCH_EPILOG, add_batch_arguments

    parser = argparse.ArgumentParser(
        description='Compute the pore geometry of a folder of structures.',
        epilog=BATCH_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('cif_folder', type=str,
                        help='folder of cif files')

    parser.add_argument('-pr', '--probe_radius', default=1.86, type=float,
                        help='probe radius (default: 1.86)')

    parser.add_argument('-ns', '--number_of_steps', default=10000, type=int,
                        help='Number of Monte Carlo samples (default: 10000)')

    parser.add_argument('-rf', '--rad_file', default=None, type=str,
                        help='path to radii file (default: None). rad file '
                             'must have .rad file extension')
    parser.add_argument('-a', '--accuracy', default='high',
                        choices=['high', 'low'],
                        help='zeo++ accuracy (default: high). low uses the '
                             'cheaper Voronoi decomposition, which is worth '
                             'trying on a structure that will not finish')
    parser.add_argument('-t', '--timeout', default=porosity_timeout,
                        type=float,
                        help='seconds allowed per structure before zeo++ is '
                             'killed and the structure recorded as a timeout '
                             f'(default: {porosity_timeout})')
    parser.add_argument('-s', '--save_dir',
                        default=read_write.DEFAULT_SAVE_DIR,
                        help='directory to save output files')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='print verbose output')
    add_batch_arguments(parser, max_time=porosity_timeout + 600)
    args = parser.parse_args()
    cif_files = [os.path.join(args.cif_folder, f) for f in os.listdir(
        args.cif_folder) if f.endswith('.cif') and not f.startswith('.')]
    compile_data(cif_files, args.save_dir, args.probe_radius,
                 args.number_of_steps, args.rad_file, args.verbose,
                 args.timeout, args.accuracy == 'high',
                 workers=args.workers, max_time=args.max_time,
                 memory_limit=args.memory_limit, shard=args.shard,
                 retry_failed=args.retry_failed, recycle=args.recycle,
                 merge_every=args.merge_every)


if __name__ == '__main__':
    main()
