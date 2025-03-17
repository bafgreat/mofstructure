#!/usr/bin/python
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
import os
import argparse
import pandas as pd
import numpy as np
import shutil
import json
from ase.io import read
import mofstructure.mofdeconstructor as MOF_deconstructor
from mofstructure.porosity import zeo_calculation
import mofstructure.filetyper as read_write


def remove_t_factor_and_unique(metal_sites):
    unique_sites = []
    seen = set()

    for site in metal_sites:
        site_without_t_factor = {key: value for key, value in site.items(
        ) if key not in ['t_factor', 'unique', 'problematic', 'type']}
        site_tuple = tuple(sorted(site_without_t_factor.items()))
        if site_tuple not in seen:
            unique_sites.append(site_without_t_factor)
            seen.add(site_tuple)

    return unique_sites


def convert_numpy_types(data):
    if isinstance(data, dict):
        return {key: convert_numpy_types(value) for key, value in data.items()}
    elif isinstance(data, list):
        return [convert_numpy_types(element) for element in data]
    elif isinstance(data, np.integer):
        return int(data)
    elif isinstance(data, np.floating):
        return float(data)
    else:
        return data


def remove_guest(ase_atom):
    '''
    Simple function to remove guest molecules in porous system
    '''
    index_non_guest = MOF_deconstructor.remove_unbound_guest(ase_atom)
    return ase_atom[index_non_guest]


def merge_two_dicts(dict_1, dict_2):
    '''
    Simple fucntion for mergeing two dictionary
    '''
    combined_dict = dict_1.copy()   # start with keys and values of x
    combined_dict.update(dict_2)    # modifies z with keys and values of y
    return combined_dict


def compile_data(cif_files, result_folder, probe_radius=1.86,\
    number_of_steps=10000, rad_file=None, verbose=False):
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
    porosity_dic = {}
    ase_atoms_dic = {}
    seen = []
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)
        # porosity_dic = {}
    else:
        try:
            porosity_dic = read_write.load_data(
                result_folder+'/porosity_dic.json')
            # porosity_dic = json.loads(pd.read_csv(result_folder+'/porosity_data.csv', index_col=False).to_json(orient='records'))
            seen = list(porosity_dic.keys())
        except Exception:
            pass

    for cif_file in cif_files:

        try:
            tmp_metal = {}
            base_name = cif_file[:cif_file.rindex('.')].split('/')[-1]
            if not base_name in seen:
                ase_atom = read(cif_file)
                ase_atom = remove_guest(ase_atom)
                if rad_file is None:
                    pores = zeo_calculation(ase_atom, probe_radius, number_of_steps)
                else:
                    pores = zeo_calculation(
                        ase_atom, probe_radius=probe_radius, number_of_steps=number_of_steps, rad_file=rad_file)
                porosity_dic[base_name] = pores
                porosity_dic = convert_numpy_types(porosity_dic)
                read_write.append_json(
                    porosity_dic, result_folder+'/porosity_data.json')

        except Exception:
            print(f"Error processing {cif_file}")
            if verbose:
                print(f"Error processing {cif_file}")
                continue
            else:
                raise Exception(f"Error processing {cif_file}")
            pass

    data_f = pd.DataFrame.from_dict(porosity_dic, orient='index')
    data_f.index.name = 'mof_names'
    data_f.to_csv(result_folder+'/porosity_data.csv')

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
