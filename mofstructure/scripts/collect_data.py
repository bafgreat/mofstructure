#!/usr/bin/python
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
import os
import argparse
import pandas as pd
from mofstructure import structure
import mofstructure.filetyper as read_write

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
    structural_data = {}
    data_to_json['sbu_smile'] = []
    data_to_json['sbu_inchikey'] = []
    data_to_json['sbu_inchi'] = []
    data_to_json['n_sbu_point_of_extension'] = []
    data_to_json['sbu_type'] = []
    data_to_json['linker_smile'] = []
    data_to_json['linker_inchikey'] = []
    data_to_json['linker_inchi'] =[]
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
    structural_data = {}
    path_to_file = f'{xyz_path}/{base_name}'
    if not os.path.exists(path_to_file):
        os.makedirs(path_to_file )
    data_to_json['n_ligands'] = len(organic_ligands)
    smi, inchikey, inchi, iupac = [], [], [], []
    for j, mof_ligand in enumerate(organic_ligands):
        smi.append(mof_ligand.info['smi'])
        inchikey.append(mof_ligand.info['inchikey'])
        iupac.append(iupacnames.get(mof_ligand.info['inchikey', None]))
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


def compile_data(cif_files, result_folder, verbose=False, oms=False):
    '''
    A workflow to remove guest, compute porosity and deconstructure
    mofs and creates a MOF database. The function starts with checking and
    removing any unbound guest molecule present in the MOF. After that it
    computed the porosity of all the MOFs and load them in a single csv.
    Finally the deconstructs the MOFs into the various building units and
    creates a three json files

    1. ase_atoms_building_units.json
    Json file containing ase atom object of all the building uints and for
    each ase atom object there are additional information in the info[] key.

    2. sbus_and_linkers.json
    A json file containing all the information about the linkers and metal
    sbu.

    3. cluster_and_ligands.json
    A json file containing all the information about the ligands and metal
    cluster.

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
    xyz_path = os.path.join(result_folder, "XYZ_DB")
    if not os.path.exists(xyz_path):
        os.makedirs(xyz_path)

    path2sbu = os.path.join(structure_db, 'sbus_and_linkers.json')
    path2ligand = os.path.join(structure_db, 'ligands_data.json')
    porosity_path =  os.path.join(structure_db, 'porosity_data.json')
    oms_path =  os.path.join(structure_db, 'structure_oms_and_general_info.json')

    if os.path.exists(path2sbu):
        all_sbu_data = read_write.load_data(path2sbu)
    else:
        all_sbu_data = {}

    if os.path.exists(path2ligand):
        all_ligand_data = read_write.load_data(path2ligand)
    else:
        all_ligand_data = {}

    if os.path.exists(porosity_path):
        all_porosity_data = read_write.load_data(porosity_path)
    else:
        all_porosity_data = {}

    if os.path.exists(oms_path):
        all_oms_data = read_write.load_data(oms_path)
    else:
        all_oms_data = {}


    seen = list(all_sbu_data.keys())

    for cif_file in cif_files:
        try:
            mof_object = structure.MOFstructure(filename=cif_file)

            base_name = os.path.basename(cif_file).split('.')[0]
            if base_name not in seen:
                print("======================================\n")
                print(f'     processing : {base_name}     \n')
                print("======================================")
                sbu_data = mof_object.get_sbu()
                ligand_data = mof_object.get_ligands()

                if sbu_data is not None:
                    metal_sbus, organic_sbus = sbu_data
                    data_to_json = collect_sbus(metal_sbus,
                                                organic_sbus,
                                                base_name,
                                                xyz_path)
                    all_sbu_data[base_name] = data_to_json
                    read_write.append_json(all_sbu_data, path2sbu)

                if ligand_data is not None:
                    metal_clusters, organic_ligands = ligand_data

                    data_to_json = collect_ligand(organic_ligands, base_name, xyz_path)
                    all_ligand_data[base_name] = data_to_json
                    read_write.append_json(all_ligand_data, path2ligand)
                porosity = mof_object.get_porosity()

                if len(porosity) > 0:
                    all_porosity_data[base_name] = porosity
                    read_write.append_json(all_porosity_data, porosity_path)
                else:
                    all_porosity_data[base_name]= None

                if len(mof_object.ase_atoms)> 5000:
                    print('system size too large, will run out of application memory')
                    print('so will skip')
                    all_oms_data[base_name] = None
                    continue
                all_oms_data[base_name] = mof_object.get_oms()
                read_write.append_json(all_oms_data, oms_path)
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
    Command line interface to deconstruct MOFs to building units,
    compute porosity and open metal sites
    '''
    parser = argparse.ArgumentParser(
        description='Run work_flow function with optional verbose output')
    parser.add_argument('cif_folder', type=str,
                        help='list of cif files. like glob')

    parser.add_argument('-o', '--oms', action='store_true',
                        help='run oms')
    parser.add_argument('-s', '--save_dir', type=str,
                        default='MOFDb', help='directory to save output files')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='print verbose output')
    args = parser.parse_args()
    cif_files = [os.path.join(args.cif_folder, f) for f in os.listdir(
        args.cif_folder) if f.endswith('.cif')]
    compile_data(cif_files, args.save_dir, args.verbose, args.oms)
