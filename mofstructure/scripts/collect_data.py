#!/usr/bin/python
'''
Command line entry point for `mofstructure_database`.

Runs the full analysis over a folder of structures: guest removal, building
unit deconstruction, porosity, ligand fingerprints, and optionally open metal
sites and topology. Results accumulate in one folder, a json file per analysis
under `Structure_Data` with a csv summary beside it, and the building units as
xyz files under `XYZ_DB`. A structure already recorded is skipped, so an
interrupted run can be resumed by repeating the command.
'''
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
import os
import argparse
from mofstructure import structure
import mofstructure.filetyper as read_write
from mofstructure.porosity import DEFAULT_TIMEOUT as default_porosity_timeout
from mofstructure.mofdeconstructor import lookup_iupac_name

iupacnames = read_write.load_iupac_names()


def collect_sbus(metal_sbus, organic_sbus, base_name, xyz_path):
    '''
    Function to compile secondary building units and region of MOFs

    Parameters
    ----------
    ase_atom : ASE atoms object
    path_to_file : result directory or folder
    '''

    path_to_file = f'{xyz_path}/{base_name}'
    if not os.path.exists(path_to_file):
        os.makedirs(path_to_file)

    data_to_json = {}
    data_to_json['sbu_smile'] = []
    data_to_json['sbu_inchikey'] = []
    data_to_json['sbu_inchi'] = []
    data_to_json['n_sbu_point_of_extension'] = []
    data_to_json['sbu_type'] = []
    data_to_json['linker_smile'] = []
    data_to_json['linker_inchikey'] = []
    data_to_json['linker_inchi'] = []
    data_to_json['n_sbu_point_of_extension'] = []
    data_to_json['n_linker_point_of_extension'] = []
    data_to_json['sbu_type'] = ''
    data_to_json['n_metal_sbus'] = len(metal_sbus)
    data_to_json['n_organic_sbus'] = len(organic_sbus)
    sbu_type, smi, inchikey, inchi, point_of_extension = [], [], [], [], []

    for i, sbu_metal in enumerate(metal_sbus):
        if len(metal_sbus) > 0:
            sbu_type.append(sbu_metal.info['sbu_type'])
            smi.append(sbu_metal.info['smi'])
            inchikey.append(sbu_metal.info['inchikey'])
            inchi.append(sbu_metal.info['inchi'])
            point_of_extension.append(
                len(sbu_metal.info['point_of_extension']))
            sbu_metal.write(f'{path_to_file}/{base_name}_metal_sbu_{i+1}.xyz')

    if len(smi) > 0:
        data_to_json['sbu_smile'] = smi

    if len(inchikey) > 0:
        data_to_json['sbu_inchikey'] = inchikey

    if len(inchi) > 0:
        data_to_json['sbu_inchi'] = inchi

    if len(point_of_extension) > 0:
        data_to_json['n_sbu_point_of_extension'] = point_of_extension

    if len(sbu_type) > 0:
        data_to_json['sbu_type'] = sbu_type

    smi, inchikey, inchi, point_of_extension = [], [], [], []
    if len(organic_sbus) > 0:
        for j, sbu_linker in enumerate(organic_sbus):
            try:
                smi.append(sbu_linker.info['smi'])
                inchikey.append(sbu_linker.info['inchikey'])
                inchi.append(sbu_linker.info['inchi'])
                point_of_extension.append(
                    len(sbu_linker.info['point_of_extension']))
                sbu_linker.write(f'{path_to_file}/{base_name}_organic_sbu_{j+1}.xyz')
            except AttributeError:
                pass

    if len(smi) > 0:
        data_to_json['linker_smile'] = smi
    if len(inchikey) > 0:
        data_to_json['linker_inchikey'] = inchikey
    if len(inchi) > 0:
        data_to_json['linker_inchi'] = inchi
    if len(point_of_extension) > 0:
        data_to_json['n_linker_point_of_extension'] = point_of_extension
    return read_write.convert_numpy_types(data_to_json)


def collect_ligand(organic_ligands, base_name, xyz_path):
    '''
    Function to compile organic ligands, metal clusters and region of MOFs

    Parameters
    ----------
    ase_atom : ASE atoms object
    path_to_file : result directory or folder
    '''
    data_to_json = {}
    path_to_file = f'{xyz_path}/{base_name}'
    if not os.path.exists(path_to_file):
        os.makedirs(path_to_file)
    data_to_json['n_ligands'] = len(organic_ligands)
    smi, inchikey, inchi, iupac = [], [], [], []
    for j, mof_ligand in enumerate(organic_ligands):
        smi.append(mof_ligand.info['smi'])
        inchikey.append(mof_ligand.info['inchikey'])
        iupac.append(lookup_iupac_name(mof_ligand.info['smi'], iupacnames))
        inchi.append(mof_ligand.info['inchi'])
        mof_ligand.write(f'{path_to_file}/{base_name}_organic_ligand_{j+1}.xyz')
    if len(smi) > 0:
        data_to_json['ligand_smile'] = smi
    else:
        data_to_json['ligand_smile'] = []
    if len(iupac) > 0:
        data_to_json['ligand_names'] = iupac
    else:
        data_to_json['ligand_names'] = []
    if len(inchikey) > 0:
        data_to_json['ligand_inchikey'] = inchikey
    else:
        data_to_json['ligand_inchikey'] = []
    if len(inchi) > 0:
        data_to_json['ligand_inchi'] = inchi
    else:
        data_to_json['ligand_inchi'] = []
    return read_write.convert_numpy_types(data_to_json)


def _legacy_state(structure_db, topology):
    '''
    What a database written before the batch runner already holds.

    **returns:**
        (complete, backfill)
            Names with every requested analysis, which are skipped, and
            names with building units but no fingerprint or net, for which
            only the missing analyses are run, as the command has always done.
    '''
    def keys(name):
        path = os.path.join(structure_db, name)
        if not os.path.exists(path) or not os.path.getsize(path):
            return set()
        return {k for k, v in read_write.load_data(path).items() if v}

    sbu = keys('sbus_and_linkers.json')
    complete = sbu & keys('fingerprint_data.json')
    if topology:
        complete &= keys('topology_data.json')
    return complete, sbu - complete


def compile_data(cif_files, result_folder, verbose=False, oms=False,
                 topology=True, topology_method='auto',
                 porosity_timeout=default_porosity_timeout,
                 probe_radius=1.86, high_accuracy=True, workers=1,
                 max_time=7200, memory_limit=None, shard=None,
                 retry_failed=False, recycle=50, topology_timeout=300,
                 merge_every=60):
    '''
    Build or extend a structure database from a list of structure files.

    For each structure the guests are removed and the building units, ligand
    cluster fingerprint, porosity, optionally the open metal sites and the
    underlying net are recorded. Every structure runs in a separate worker
    process with a hard time limit (see `mofstructure.batch`), results are
    checkpointed as they finish, and the json and csv files under
    `Structure_Data` are rebuilt at the end. Repeating the call continues an
    interrupted run.

    **parameters:**
        - cif_files: list of str
            CIF files, or anything ASE reads, each holding one framework.

        - result_folder: str
            Directory the database is written to.

        - verbose: bool
            Print each structure as it is processed.

        - oms: bool
            Record open metal sites.

        - topology: bool
            Compute the underlying net.

        - topology_method: str
            Node definition passed to `MOFstructure.get_topology`, or "auto".

        - porosity_timeout: float
            Seconds allowed for zeo++ per structure.

        - probe_radius: float
            Probe radius used for the porosity.

        - high_accuracy: bool
            Run zeo++ at high accuracy.

        - workers: int
            Worker processes.

        - max_time: float
            Hard wall-clock limit per structure in seconds.

        - memory_limit: float or None
            Memory limit per worker in GB (Linux).

        - shard: str or None
            "i/N" to process one part of the list.

        - retry_failed: bool
            Run again structures recorded as failed.

        - recycle: int
            Structures per worker before it is replaced.

        - topology_timeout: int
            Seconds allowed for net identification.

        - merge_every: float
            Minutes between refreshes of the json and csv files during the
            run; 0 refreshes only at the end.
    '''
    from types import SimpleNamespace
    from mofstructure.batch import progress_printer
    from mofstructure.batch_tasks import database_name, run_command

    structure_db = os.path.join(result_folder, read_write.STRUCTURE_DATA)
    xyz_path = os.path.join(result_folder, read_write.XYZ_DB)
    os.makedirs(structure_db, exist_ok=True)
    os.makedirs(xyz_path, exist_ok=True)

    items = [(database_name(f), f) for f in cif_files]
    options = SimpleNamespace(workers=workers, max_time=max_time,
                              memory_limit=memory_limit, shard=shard,
                              retry_failed=retry_failed, recycle=recycle,
                              merge_every=merge_every)
    kwargs = dict(xyz_path=xyz_path, oms=oms, topology=topology,
                  topology_method=topology_method,
                  topology_timeout=topology_timeout,
                  porosity_timeout=porosity_timeout,
                  probe_radius=probe_radius, high_accuracy=high_accuracy)

    printer = progress_printer()

    def quiet(record, done, total):
        # Failures always, otherwise every hundredth structure, so a long
        # run leaves a readable log.
        if record['status'] != 'ok' or done == total or done % 100 == 0:
            printer(record, done, total)

    complete, backfill = _legacy_state(structure_db, topology)
    if backfill:
        # Building units recorded by an earlier version: add only the
        # fingerprint and the net, keeping the stored porosity and units.
        only = ['fingerprint'] + (['topology'] if topology else [])
        run_command('database', [i for i in items if i[0] in backfill],
                    'mofstructure.batch_tasks:database_task',
                    {**kwargs, 'only': only}, result_folder, options,
                    on_record=None if verbose else quiet,
                    extra_skip=complete)
    run_command('database', [i for i in items if i[0] not in backfill],
                'mofstructure.batch_tasks:database_task',
                kwargs, result_folder, options,
                on_record=None if verbose else quiet,
                extra_skip=complete)
    if verbose:
        print(f"Saved results to {result_folder}")


def main():
    '''
    Command line interface to deconstruct MOFs to building units,
    compute porosity, topology, fingerprints and open metal sites.
    '''
    from mofstructure.batch_tasks import BATCH_EPILOG, add_batch_arguments

    parser = argparse.ArgumentParser(
        description='Build a structure database from a folder of '
                    'structures.',
        epilog=BATCH_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('cif_folder', type=str,
                        help='folder of cif files')

    parser.add_argument('--oms', action='store_true',
                        help='run oms')
    parser.add_argument('-t', '--topology', action='store_true',
                        default=True,
                        help='compute the underlying net (on by default)')
    parser.add_argument('--no-topology', action='store_false',
                        dest='topology',
                        help='skip the net, leaving the rest of the record')
    parser.add_argument('--method', '--topology_method', type=str,
                        dest='topology_method',
                        default='auto',
                        choices=['auto', 'sbus', 'ligand_cluster',
                                 'all_node', 'single_node'],
                        help='node definition used to build the net; "auto" '
                             'chooses it from the structure, so a folder of '
                             'MOFs, COFs and zeolites is handled in one go '
                             '(--topology_method is a deprecated alias)')
    parser.add_argument('--topology_timeout', type=int, default=300,
                        help='seconds allowed for net identification '
                             '(default: 300)')
    parser.add_argument('-pr', '--probe_radius', default=1.86, type=float,
                        help='probe radius used for the porosity '
                             '(default: 1.86)')
    parser.add_argument('-a', '--accuracy', default='high',
                        choices=['high', 'low'],
                        help='zeo++ accuracy (default: high)')
    parser.add_argument('--porosity_timeout', type=float,
                        default=default_porosity_timeout,
                        help='seconds allowed per structure before zeo++ is '
                             'killed and the structure recorded as a timeout '
                             f'(default: {default_porosity_timeout})')
    parser.add_argument('-s', '--save_dir', type=str,
                        default=read_write.DEFAULT_SAVE_DIR,
                        help='directory to save output files')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='print a line for every structure')
    add_batch_arguments(parser, max_time=7200)
    args = parser.parse_args()
    cif_files = [os.path.join(args.cif_folder, f) for f in os.listdir(
        args.cif_folder) if f.endswith('.cif') and not f.startswith('.')]
    compile_data(cif_files, args.save_dir, args.verbose, args.oms,
                 args.topology, args.topology_method, args.porosity_timeout,
                 args.probe_radius, args.accuracy == 'high',
                 workers=args.workers, max_time=args.max_time,
                 memory_limit=args.memory_limit, shard=args.shard,
                 retry_failed=args.retry_failed, recycle=args.recycle,
                 merge_every=args.merge_every,
                 topology_timeout=args.topology_timeout)


if __name__ == '__main__':
    main()
