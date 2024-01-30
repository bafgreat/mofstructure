#!/usr/bin/python
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
import argparse
import os
from ase.io import read
import pandas as pd
import mofstructure.mofdeconstructor as MOF_deconstructor
from mofstructure.porosity import zeo_calculation


def sbu_data(ase_atom, path_to_file):
    '''
    Function to compile secondary building units and region of MOFs

    Parameters
    ----------
    ase_atom : ASE atoms object
    path_to_file : result directory or folder
    '''
    connected_components, atoms_indices_at_breaking_point, porpyrin_checker, all_regions = MOF_deconstructor.\
        secondary_building_units(ase_atom)
    metal_sbus, organic_sbus, building_unit_regions = MOF_deconstructor.\
        find_unique_building_units(
            connected_components,
            atoms_indices_at_breaking_point,
            ase_atom, porpyrin_checker,
            all_regions,
            cheminfo=True
        )
    for i, sbu_metal in enumerate(metal_sbus):
        # print('sbu', sbu_metal.info['sbu_type'])
        # print('point_of_extension', len(sbu_metal.info['point_of_extension']))
        sbu_metal.write(path_to_file + '_metal_sbu_'+str(i+1)+'.xyz')
    for j, sbu_linker in enumerate(organic_sbus):
        sbu_linker.write(path_to_file + '_organic_sbu_'+str(j+1)+'.xyz')
    return building_unit_regions


def ligand_data(ase_atom, path_to_file):
    '''
    Function to compile organic ligands, metal clusters and region of MOFs

    Parameters
    ----------
    ase_atom : ASE atoms object
    path_to_file : result directory or folder
    '''
    connected_components, atoms_indices_at_breaking_point, porpyrin_checker, all_regions = MOF_deconstructor.\
        ligands_and_metal_clusters(ase_atom)
    metal_clusters, organic_ligands, building_unit_regions = MOF_deconstructor.\
        find_unique_building_units(
            connected_components,
            atoms_indices_at_breaking_point,
            ase_atom, porpyrin_checker,
            all_regions,
            cheminfo=True
        )
    for i, mof_metal in enumerate(metal_clusters):
        mof_metal.write(path_to_file+'_metal_cluster_'+str(i+1)+'.xyz')
    for j,  mof_ligand in enumerate(organic_ligands):
        mof_ligand.write(path_to_file + '_organic_ligand_'+str(j+1)+'.xyz')
    return building_unit_regions


def remove_guest(ase_atom):
    '''
    Simple function to remove guest molecules in porous system
    '''
    index_non_guest = MOF_deconstructor.remove_unbound_guest(ase_atom)
    return ase_atom[index_non_guest]


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
    try:
        ase_atom = read(cif_file)
        ase_atom = remove_guest(ase_atom)
        base_name = cif_file[:cif_file.rindex('.')].split('/')[-1]
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        path_to_file = save_dir+'/'+base_name

        ase_atom = remove_guest(ase_atom)
        pores = zeo_calculation(ase_atom)
        data_f = pd.DataFrame(pores, index=[0])
        data_f.to_csv(path_to_file+'_porosity_data.csv')
        sbu_data(ase_atom, path_to_file)
        ligand_data(ase_atom, path_to_file)
    except Exception:
        pass
    if verbose:
        print(f"Saved results to {path_to_file}")
    return

def main():
    '''
    Routine runner for CLI
    '''
    parser = argparse.ArgumentParser(description='Run work_flow function with optional verbose output')
    parser.add_argument('cif_file', type=str, help='path to CIF file')
    parser.add_argument('-s', '--save_dir', type=str, default='MOF_building_units', help='directory to save output files')
    parser.add_argument('-v', '--verbose', action='store_true', help='print verbose output')
    args = parser.parse_args()
    work_flow(args.cif_file, args.save_dir, args.verbose)
