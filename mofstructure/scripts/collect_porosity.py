#!/usr/bin/python
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
import os
import argparse
import pandas as pd
import json
from mofstructure import structure
import mofstructure.filetyper as read_write



def compile_data(cif_files,
                 result_folder,
                 probe_radius=1.86,
                 number_of_steps=10000,
                 rad_file=None,
                 verbose=False):
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
    structure_db = os.path.join(result_folder, "Structure_Data")
    if not os.path.exists(structure_db):
        os.makedirs(structure_db)
    porosity_path =  os.path.join(structure_db, 'porosity_data.json')
    if os.path.exists(porosity_path):
        all_porosity_data = read_write.load_data(porosity_path)
    else:
        all_porosity_data = {}

    seen = list(all_porosity_data.keys())

    for cif_file in cif_files:
        try:
            mof_object = structure.MOFstructure(filename=cif_file)

            base_name = os.path.basename(cif_file).split('.')[0]
            if base_name not in seen:
                print("======================================\n")
                print(f'     processing : {base_name}     \n')
                print("======================================")

                porosity = mof_object.get_porosity()

                if len(porosity) > 0:
                    all_porosity_data[base_name] = porosity
                    read_write.append_json(all_porosity_data, porosity_path)
                else:
                    all_porosity_data[base_name]= None
            else:
                print("======================================\n")
                print(f" !!! {base_name} is already done !!! \n")
                print("======================================")
        except Exception:
            pass

    data_f = pd.DataFrame.from_dict(all_porosity_data, orient='index')
    data_f.index.name = 'mof_names'
    data_f.to_csv(structure_db+'/porosity_data.csv')

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
    parser.add_argument('-s', '--save_dir',
                        default='MOFDb', help='directory to save output files')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='print verbose output')
    args = parser.parse_args()
    cif_files = [os.path.join(args.cif_folder, f) for f in os.listdir(
        args.cif_folder) if f.endswith('.cif')]
    compile_data(cif_files, args.save_dir, args.probe_radius,
                 args.number_of_steps,  args.rad_file, args.verbose)
