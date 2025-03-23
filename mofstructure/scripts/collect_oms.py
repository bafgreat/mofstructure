#!/usr/bin/python
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
import os
import argparse
import pandas as pd
from mofstructure import filetyper
from mofstructure import structure
import mofstructure.filetyper as read_write


def compile_data(cif_files, result_folder, verbose=False):
    '''
    A workflow to remove guest and compute open metal sites from a folder containing
    cif files and create a database of open metal sites. This code can work for all
    periodic systems containing metals.

    The function starts with checking and removing any unbound
    guest molecule present in the periodic system.

    1. metal_info.json
    A json file containing all information about the metals in the system.
    This file can easily be converted to csv format.

    ::
        Parameters
        ----------
        cif_file : a cif file or any ase readable file containing a MOF.
        result_folder : path to output folder
    '''
    metal_info = {}
    seen = []
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)
    else:
        try:
            metal_info = read_write.load_data(result_folder+'/structure_oms_and_general_info.json')
            seen = list(metal_info.keys())
        except Exception:
            pass


    for cif_file in cif_files:
        try:
            base_name = os.path.basename(cif_file).split('.')[0]
            print("======================================\n")
            print(f'     processing : {base_name}     \n')
            print("======================================")
            if base_name not in seen:
                print(f'processing: {base_name}')
                mof_object = structure.MOFstructure(filename=cif_file)
                if len(mof_object.ase_atoms)> 5000:
                    print('system size too large, will run out of application memory')
                    print('so will skip')
                    continue
                oms_data = mof_object.get_oms()

                metal_info[base_name] = oms_data
                read_write.append_json(metal_info, result_folder+'/structure_oms_and_general_info.json')
            else:
                print("======================================\n")
                print(f" Already Done {base_name}\n")
                print("======================================")
        except Exception:
            print(f'failed processing {base_name}')
            pass

    data_f = pd.DataFrame.from_dict(metal_info, orient='index')
    data_f.to_csv(result_folder+'/structure_oms_and_general_info.csv')

    if verbose:
        print(f"Saved results to {result_folder}")
    return


def main():
    '''
    mofstructure command line interface for computiong open metal sites
    '''
    parser = argparse.ArgumentParser(
        description='Run work_flow function with optional verbose output')
    parser.add_argument('cif_folder', type=str,
                        help='list of cif files. like glob')


    parser.add_argument('-s', '--save_dir', type=str,
                        default='MOFDb', help='directory to save output files')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='print verbose output')
    args = parser.parse_args()
    if os.path.isdir(args.cif_folder):
        cif_files = [os.path.join(args.cif_folder, f)
                 for f in os.listdir(args.cif_folder)
                 if f.endswith('.cif')]
    else:
        cif_files = [args.cif_folder]
    compile_data(cif_files, args.save_dir, args.verbose)
