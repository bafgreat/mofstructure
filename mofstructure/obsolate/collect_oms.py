#!/usr/bin/python
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
import os
import argparse
import pandas as pd
import shutil
from ase.io import read
from omsdetector_forked import MofCollection
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


def find_oms(cif_file, base_name, ase_atom, result_folder):
    """
    Function to compute open metal sites
    """
    result = []
    oms = []
    enviroment = {}
    test_folder = result_folder+'/test'
    a_mof_collection = MofCollection(
        path_list=[cif_file], analysis_folder=test_folder)
    print ("This is important")
    a_mof_collection.analyse_mofs(overwrite=True)
    data = read_write.load_data(
        test_folder+'/oms_results/' + base_name+'/'+base_name+'.json')

    # result['has_oms'] = data['has_oms']
    metal_environment = MOF_deconstructor.metal_coordination_enviroment(
        ase_atom)
    metal_sites = data['metal_sites']
    unique_metal_sites = remove_t_factor_and_unique(metal_sites)

    # cn = {metal_site["metal"]:metal_site["number_of_linkers"] for metal_site in metal_sites}
    # print ( "coordination number", cn)
    # max_cn = max(list(cn.values()))
    # print (max_cn , cn )
    # result['metal_cn'] = cn
    # result['max_cn'] = max_cn
    # if data['has_oms'] is True:

    for dat in unique_metal_sites:
        # if dat['is_open'] is True:
        # result["metal"] = dat['metal']
        tmp = {}
        tmp["metal"] = dat['metal']
        # tmp["type"] = dat["type"]
        tmp["is_open"] = dat["is_open"]
        tmp['cn'] = dat["number_of_linkers"]
        tmp['environment'] = metal_environment[dat['metal']]
        result.append(tmp)

        # oms.append(metal)
        # enviroment[metal] = metal_environment[metal]
    # result['oms'] = list(set(oms))
    # result['environment'] = enviroment
    # else:
    #     result['oms'] = []
    #     result['environment'] = {}
    # if os.path.exists(test_folder):
    #     shutil.rmtree(test_folder)
    return result


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
        # porosity_dic = {}
    else:
        try:
            metal_info = read_write.load_data(result_folder+'/metal_info.json')
            seen = list(metal_info.keys())
        except Exception:
            pass

    for cif_file in cif_files:
        # try:
        tmp_metal = {}
        # encoder = read_write.AtomsEncoder
        base_name = cif_file[:cif_file.rindex('.')].split('/')[-1]
        ase_atom = read(cif_file)
        ase_atom = remove_guest(ase_atom)
        if not base_name in seen:

            open_metal_sites = find_oms(cif_file, base_name,  ase_atom, result_folder)
            print (open_metal_sites)
            metal_info[base_name] = open_metal_sites
            read_write.append_json(metal_info, result_folder+'/metal_info.json')
        # except Exception:
        #     pass

    data_f = pd.DataFrame.from_dict(metal_info, orient='index')
    # data_f.index.name = 'mof_names'
    data_f.to_csv(result_folder+'/metal_info.csv')

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
    cif_files = [os.path.join(args.cif_folder, f) for f in os.listdir(
        args.cif_folder) if f.endswith('.cif')]
    compile_data(cif_files, args.save_dir, args.verbose)
