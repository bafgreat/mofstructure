#!/usr/bin/python
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
import argparse
import os
import pandas as pd
from mofstructure import filetyper
from mofstructure import structure


def collect_sbus(metal_sbus, organic_sbus, base_name, xyz_path):

    for i, sbu_metal in enumerate(metal_sbus):
        sbu_metal.write(f'{xyz_path}/{base_name}_metal_sbu_{i+1}.xyz')
    for j, sbu_linker in enumerate(organic_sbus):
        sbu_linker.write(f'{xyz_path}/{base_name}_organic_sbu_{j+1}.xyz')


def collect_ligand(organic_ligands, base_name, xyz_path):

    for j, sbu_linker in enumerate(organic_ligands):
        sbu_linker.write(f'{xyz_path}/{base_name}_organic_ligand_{j+1}.xyz')


def work_flow(cif_file, save_dir, verbose=False):
    '''
    A workflow to remove guest, compute porosity and deconstructure
    mofs. The function starts with checking and removing any unbound
    guest molecule present in the MOF. After that it computed the porosity
    of the MOF and load this data in the output folder. Finally the deconstructs
    the MOFs into the various building units and safe these building units
    in '.xyz' formats in the output folder.
    Parameters
    ----------
    cif_file : a cif file or any ase readable file containing a MOF.
    save_dir : path to output folder
    '''
    # try:
    mof_object = structure.MOFstructure(filename=cif_file)

        # base_name = cif_file[:cif_file.rindex('.')].split('/')[-1]

    base_name = os.path.basename(cif_file).split('.')[0]
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    path_to_file = save_dir+'/'+base_name

    xyz_path = save_dir+'/XYZ_DB'
    if not os.path.exists(xyz_path):
        os.makedirs(xyz_path)

    guest_free_atom = mof_object.remove_guest()
    guest_free_atom.write(f'{path_to_file}_guest_free.cif')

    sbu_data = mof_object.get_sbu()
    if sbu_data is not None:
        metal_sbus, organic_sbus = sbu_data
        collect_sbus(metal_sbus, organic_sbus, base_name, xyz_path)

    cluster_ligand = mof_object.get_ligands()
    if cluster_ligand is not None:
        _, organic_ligands = cluster_ligand
        collect_ligand(organic_ligands, base_name, xyz_path)

    pores = mof_object.get_porosity()
    filetyper.write_json(pores, path_to_file+'_porosity_data.json')

    data_f = pd.DataFrame(pores, index=[0])
    data_f.to_csv(path_to_file+'_porosity_data.csv')

    oms_data = mof_object.get_oms()
    filetyper.write_json(oms_data, path_to_file+'_general_info.json')

    # except Exception:
    #     pass
    if verbose:
        print(f"Saved results to {path_to_file}")
    return


def main():
    '''
    Routine runner for CLI
    '''
    parser = argparse.ArgumentParser(description='Run work_flow function with optional verbose output')
    parser.add_argument('cif_file', type=str, help='path to CIF file')
    parser.add_argument('-s', '--save_dir', type=str, default='mof_building_units', help='directory to save output files')
    parser.add_argument('-v', '--verbose', action='store_true', help='print verbose output')
    args = parser.parse_args()
    work_flow(args.cif_file, args.save_dir, args.verbose)
