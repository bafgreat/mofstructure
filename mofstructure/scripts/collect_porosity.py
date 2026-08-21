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
from mofstructure import structure
import mofstructure.filetyper as read_write
from mofstructure.porosity import DEFAULT_TIMEOUT as porosity_timeout
from mofstructure.porosity import empty_porosity_record


def compile_data(cif_files,
                 result_folder,
                 probe_radius=1.86,
                 number_of_steps=10000,
                 rad_file=None,
                 verbose=False,
                 timeout=porosity_timeout,
                 high_accuracy=True):
    '''
    A workflow to remove guest and compute porosity from any porous periodic system.
    The results is written in both a json format and csv file format.
    The function starts with checking and removing any unbound
    guest molecule present in the porous. After that it computed the porosity
    of all the systems and load them in a single csv. The function always computes
    the high accuracy calculation.

    1. ase_atoms_building_units.json
    Json file containing ase atom object of all the building uints and for
    each ase atom object there are additional information in the info[] key.

    2. sbus_and_linkers.json
    A json file containing all the information about the linkers and metal
    sbu.

    3. cluster_and_ligands.json
    A json file containing all the information about the ligands and metal
    cluster.
    ::
        Parameters
        ----------
        cif_file : a cif file or any ase readable file containing a MOF.
        result_folder : path to output folder
    '''
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)
    structure_db = os.path.join(result_folder, read_write.STRUCTURE_DATA)
    if not os.path.exists(structure_db):
        os.makedirs(structure_db)
    porosity_path = os.path.join(structure_db, 'porosity_data.json')
    if os.path.exists(porosity_path):
        all_porosity_data = read_write.load_data(porosity_path)
    else:
        all_porosity_data = {}

    seen = list(all_porosity_data.keys())

    for cif_file in cif_files:
        base_name = os.path.basename(cif_file).split('.')[0]
        try:
            mof_object = structure.MOFstructure(filename=cif_file)
            if base_name not in seen:
                print("======================================\n")
                print(f'     processing : {base_name}     \n')
                print("======================================")

                porosity = mof_object.get_porosity(
                    probe_radius=probe_radius,
                    number_of_steps=number_of_steps,
                    rad_file=rad_file,
                    timeout=timeout,
                    high_accuracy=high_accuracy)
                # The record has the same keys whether or not zeo++ managed
                # the structure, so a failure is a row of missing values
                # carrying its reason rather than a gap in the table.
                all_porosity_data[base_name] = porosity
                read_write.append_json(all_porosity_data, porosity_path)
            else:
                print("======================================\n")
                print(f" !!! {base_name} is already done !!! \n")
                print("======================================")
        except Exception as error:
            # A structure that cannot even be read still belongs in the
            # table, with the reason attached. Skipping it silently leaves
            # the caller unable to tell a missing row from a failed one.
            print(f'!!! {base_name} could not be analysed: '
                  f'{type(error).__name__}: {error} !!!')
            all_porosity_data[base_name] = empty_porosity_record(
                f'failed:{type(error).__name__}')

    read_write.summary_frame(all_porosity_data).to_csv(
        structure_db + '/porosity_data.csv')

    if verbose:
        print(f"Saved results to {result_folder}")
    return


def main():
    '''
    mofstructure command line interface to compute the porosity of any periodic
    system.
    '''
    parser = argparse.ArgumentParser(
        description='Run work_flow function with optional verbose output')
    parser.add_argument('cif_folder', type=str,
                        help='list of cif files. like glob')

    parser.add_argument('-pr', '--probe_radius', default=1.86, type=float,
                        help='probe radius (default: 1.86)')

    parser.add_argument('-ns', '--number_of_steps', default=10000, type=int,
                        help='Number of GCMC simulation cycles (default: 10000)')

    parser.add_argument('-rf', '--rad_file', default=None, type=str,
                        help='path to radii file (default: None). rad file must have .rad file extension')
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
    args = parser.parse_args()
    cif_files = [os.path.join(args.cif_folder, f) for f in os.listdir(
        args.cif_folder) if f.endswith('.cif')]
    compile_data(cif_files, args.save_dir, args.probe_radius,
                 args.number_of_steps, args.rad_file, args.verbose,
                 args.timeout, args.accuracy == 'high')
