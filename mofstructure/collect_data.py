#!/usr/bin/python
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__version__ = '0.1.0'
__status__ = "production"
import sys
import os
import glob
from ase.io import read
import mofstructure.mofdeconstructor as MOF_deconstructor
from mofstructure.porosity import zeo_calculation
import mofstructure.filetyper as Writer
import argparse
import pandas as pd

def sbu_data(ase_atom):
    '''
    Function to compile secondary building units and region of MOFs

    Parameters
    ----------
    ase_atom : ASE atoms object
    path_to_file : result directory or folder
    '''
    
    data_to_json = {}
    sbu_data = {}
    connected_components, atoms_indices_at_breaking_point, porpyrin_checker, all_regions = MOF_deconstructor.\
        secondary_building_units(ase_atom)
        
    if len(connected_components) > 1:
        metal_sbus, organic_sbus, building_unit_regions = MOF_deconstructor.\
            find_unique_building_units(
                connected_components,
                atoms_indices_at_breaking_point,
                ase_atom, porpyrin_checker,
                all_regions,
                cheminfo=True
            )
            
        metal_elt, metal_coordination = MOF_deconstructor.metal_coordination_number(ase_atom)
        
        data_to_json['n_linkers'] = len(organic_sbus)
        data_to_json['n_sbus'] = len(metal_sbus)
        data_to_json['n_metals'] = len(metal_elt)
        data_to_json['max_cn'] = max(list(metal_coordination.values()))
        
        data_to_json['metals'] = list(metal_elt)
        data_to_json['metal_cn'] = metal_coordination
        sbu_type, smi, inchikey, inchi, point_of_extension  = [], [], [], [], []
        for i, sbu_metal in enumerate(metal_sbus):
            sbu_type.append(sbu_metal.info['sbu_type'])
            smi.append(sbu_metal.info['smi'])
            inchikey.append(sbu_metal.info['inchikey'])
            inchi.append(sbu_metal.info['inchi'])
            point_of_extension.append(len(sbu_metal.info['point_of_extension']))
            sbu_data[sbu_metal.info['inchikey']] = sbu_metal
            
        data_to_json['sbu_smile'] = smi
        data_to_json['sbu_inchikey'] = inchikey
        data_to_json['sbu_inchi'] = inchi
        data_to_json['n_sbu_point_of_extension'] = point_of_extension
        data_to_json['sbu_type'] = sbu_type
        
        smi, inchikey, inchi, point_of_extension  = [], [], [], []
        for j, sbu_linker in enumerate(organic_sbus):
            smi.append(sbu_linker.info['smi'])
            inchikey.append(sbu_linker.info['inchikey'])
            inchi.append(sbu_linker.info['inchi'])
            point_of_extension.append(len(sbu_linker.info['point_of_extension']))
            sbu_data[sbu_linker.info['inchikey']] = sbu_linker
            
        data_to_json['linker_smile'] = smi
        data_to_json['linker_inchikey'] = inchikey
        data_to_json['linker_inchi'] = inchi
        data_to_json['n_linker_point_of_extension'] = point_of_extension
            
    return data_to_json,  sbu_data


def ligand_data(ase_atom):
    '''
    Function to compile organic ligands, metal clusters and region of MOFs

    Parameters
    ----------
    ase_atom : ASE atoms object
    path_to_file : result directory or folder
    '''
    
    data_to_json = {}
    sbu_data = {}
    connected_components, atoms_indices_at_breaking_point, porpyrin_checker, all_regions = MOF_deconstructor.\
        ligands_and_metal_clusters(ase_atom)
    if len(connected_components) > 1:
        metal_clusters, organic_ligands, building_unit_regions = MOF_deconstructor.\
            find_unique_building_units(
                connected_components,
                atoms_indices_at_breaking_point,
                ase_atom, porpyrin_checker,
                all_regions,
                cheminfo=True
            )

        metal_elt, metal_coordination = MOF_deconstructor.metal_coordination_number(ase_atom)
        
        data_to_json['n_ligands'] = len(organic_ligands)
        data_to_json['n_clusters'] = len(metal_clusters)
        data_to_json['n_metals'] = len(metal_elt)
        data_to_json['max_cn'] = max(list(metal_coordination.values()))
        
        data_to_json['metals'] = metal_elt
        data_to_json['metal_cn'] = metal_coordination
        
        smi, inchikey, inchi = [], [], []
        for i, mof_metal in enumerate(metal_clusters):
            smi.append(mof_metal.info['smi'])
            inchikey.append(mof_metal.info['inchikey'])
            inchi.append(mof_metal.info['inchi'])
            sbu_data[mof_metal.info['inchikey']] = mof_metal
            
        data_to_json['metal_cluster_smile'] = smi
        data_to_json['metal_cluster_inchikey'] = inchikey
        data_to_json['metal_cluster_inchi'] = inchi
        
            
        smi, inchikey, inchi = [], [], []
        for j,  mof_ligand in enumerate(organic_ligands):
            smi.append(mof_ligand.info['smi'])
            inchikey.append(mof_ligand.info['inchikey'])
            inchi.append(mof_ligand.info['inchi'])
            sbu_data[mof_ligand.info['inchikey']] = mof_ligand
        
        data_to_json['ligand_smile'] = smi
        data_to_json['ligand_inchikey'] = inchikey
        data_to_json['ligand_inchi'] = inchi
    return data_to_json, sbu_data


def remove_guest(ase_atom):
    '''
    Simple function to remove guest molecules in porous system 
    '''
    index_non_guest = MOF_deconstructor.remove_unbound_guest(ase_atom)
    return ase_atom[index_non_guest]

def merge_two_dicts(x, y):
    z = x.copy()   # start with keys and values of x
    z.update(y)    # modifies z with keys and values of y
    return z
    
def compile_data(cif_files, result_folder, verbose=False):
    porosity_dic = {}
    ase_dic = {}
    search_data1 = {}
    search_data2 = {}
    
    for cif_file in cif_files:
        ase_atom = read(cif_file)
        ase_atom = remove_guest(ase_atom)
        base_name = cif_file[:cif_file.rindex('.')].split('/')[-1]
        print (base_name)
        ase_atom = remove_guest(ase_atom)
        pores = zeo_calculation(ase_atom)
       
        data_to_json1, ase_data1 = sbu_data(ase_atom)
        data_to_json2, ase_data2  = ligand_data(ase_atom)
        porosity_dic[base_name] = pores
        search_data1[base_name] = data_to_json1
        search_data2[base_name] = data_to_json2
        sbu_data3 = dict(merge_two_dicts(ase_data1,ase_data2))
        ase_dic[base_name]  = sbu_data3
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)
        
    df = pd.DataFrame.from_dict(porosity_dic, orient='index')
    df.to_csv(result_folder+'/porosity_data.csv')
    
    encoder = Writer.AtomsEncoder
    Writer.Write_JSON_ATOM(ase_dic, encoder, result_folder+'/ase_atoms_building_units.json')
    Writer.Write_Json(search_data1,  result_folder+'/sbus_and_linkers.json')
    Writer.Write_Json(search_data2, result_folder+'/cluster_and_ligands.json')
    if verbose:
        print(f"Saved results to {path_to_file}")
    return
    
def main():
    parser = argparse.ArgumentParser(description='Run work_flow function with optional verbose output')
    parser.add_argument('cif_folder', type=str, help='list of cif files. like glob')
    parser.add_argument('-s', '--save_dir', type=str, default='MOFDb', help='directory to save output files')
    parser.add_argument('-v', '--verbose', action='store_true', help='print verbose output')
    args = parser.parse_args()
    cif_files = [os.path.join(args.cif_folder, f) for f in os.listdir(args.cif_folder) if f.endswith('.cif')]
    print (cif_files)
    compile_data(cif_files, args.save_dir, args.verbose)


