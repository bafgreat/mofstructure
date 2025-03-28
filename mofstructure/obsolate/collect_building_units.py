#!/usr/bin/python
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
import os
import argparse
from ase.io import read
import mofstructure.mofdeconstructor as MOF_deconstructor
import mofstructure.filetyper as read_write


def sbu_data(ase_atom):
    '''
    Function to compile secondary building units and region of MOFs

    Parameters
    ----------
    ase_atom : ASE atoms object
    path_to_file : result directory or folder
    '''
    data_to_json = {}
    structural_data = {}
    connected_components, atoms_indices_at_breaking_point, porpyrin_checker, all_regions = MOF_deconstructor.\
        secondary_building_units(ase_atom)

    if len(connected_components) > 1:
        metal_sbus, organic_sbus, _ = MOF_deconstructor.\
            find_unique_building_units(
                connected_components,
                atoms_indices_at_breaking_point,
                ase_atom, porpyrin_checker,
                all_regions,
                cheminfo=True,
                add_dummy=True
            )
        metal_elt, metal_coordination = MOF_deconstructor.metal_coordination_number(
            ase_atom)
        data_to_json['n_linkers'] = len(organic_sbus)
        data_to_json['n_sbus'] = len(metal_sbus)
        data_to_json['n_metals'] = len(metal_elt)
        data_to_json['max_cn'] = max(list(metal_coordination.values()))
        data_to_json['metals'] = list(metal_elt)
        data_to_json['metal_cn'] = metal_coordination
        sbu_type, smi, inchikey, inchi, point_of_extension = [], [], [], [], []
        for sbu_metal in metal_sbus:
            # mof_structure  = MofStructure(lattice=sbu_metal.get_cell().tolist(), species=sbu_metal.get_chemical_symbols(), coords=sbu_metal.get_positions(), coords_are_cartesian=True)
            # # output_folder = 'check'
            # # print ('Summary ************')
            # print (MofStructure)
            # mof_structure.analyze_metals2()
            # print ('Summary ************')
            # print(mof_structure.summary)
            # print ('Open metal', metal_site.is_open())
            sbu_type.append(sbu_metal.info['sbu_type'])
            smi.append(sbu_metal.info['smi'])
            inchikey.append(sbu_metal.info['inchikey'])
            inchi.append(sbu_metal.info['inchi'])
            point_of_extension.append(
                len(sbu_metal.info['point_of_extension']))
            structural_data[sbu_metal.info['inchikey']] = sbu_metal
        data_to_json['sbu_smile'] = smi
        data_to_json['sbu_inchikey'] = inchikey
        data_to_json['sbu_inchi'] = inchi
        data_to_json['n_sbu_point_of_extension'] = point_of_extension
        data_to_json['sbu_type'] = sbu_type
        smi, inchikey, inchi, point_of_extension = [], [], [], []
        for sbu_linker in organic_sbus:
            smi.append(sbu_linker.info['smi'])
            inchikey.append(sbu_linker.info['inchikey'])
            inchi.append(sbu_linker.info['inchi'])
            point_of_extension.append(
                len(sbu_linker.info['point_of_extension']))
            structural_data[sbu_linker.info['inchikey']] = sbu_linker
        data_to_json['linker_smile'] = smi
        data_to_json['linker_inchikey'] = inchikey
        data_to_json['linker_inchi'] = inchi
        data_to_json['n_linker_point_of_extension'] = point_of_extension
    return data_to_json, structural_data


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


def ligand_data(ase_atom):
    '''
    Function to compile organic ligands, metal clusters and region of MOFs

    Parameters
    ----------
    ase_atom : ASE atoms object
    path_to_file : result directory or folder
    '''

    data_to_json = {}
    structural_data = {}
    connected_components, atoms_indices_at_breaking_point, porpyrin_checker, all_regions = MOF_deconstructor.\
        ligands_and_metal_clusters(ase_atom)
    if len(connected_components) > 1:
        metal_clusters, organic_ligands, _ = MOF_deconstructor.\
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
        # data_to_json['max_cn'] = max(list(metal_coordination.values()))
        data_to_json['metals'] = metal_elt
        # data_to_json['metal_cn'] = metal_coordination
        smi, inchikey, inchi = [], [], []
        for mof_metal in metal_clusters:
            smi.append(mof_metal.info['smi'])
            inchikey.append(mof_metal.info['inchikey'])
            inchi.append(mof_metal.info['inchi'])
            structural_data[mof_metal.info['inchikey']] = mof_metal
        data_to_json['metal_cluster_smile'] = smi
        data_to_json['metal_cluster_inchikey'] = inchikey
        data_to_json['metal_cluster_inchi'] = inchi
        smi, inchikey, inchi = [], [], []
        for mof_ligand in organic_ligands:
            smi.append(mof_ligand.info['smi'])
            inchikey.append(mof_ligand.info['inchikey'])
            inchi.append(mof_ligand.info['inchi'])
            structural_data[mof_ligand.info['inchikey']] = mof_ligand
        data_to_json['ligand_smile'] = smi
        data_to_json['ligand_inchikey'] = inchikey
        data_to_json['ligand_inchi'] = inchi
    return data_to_json, structural_data


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


def compile_data(cif_files, result_folder, verbose=False):
    '''
    A workflow to remove guest, compute porosity and deconstructure
    mofs and creates a MOF database. The function starts with checking and removing any unbound
    guest molecule present in the MOF. After that it computed the porosity
    of all the MOFs and load them in a single csv. Finally the deconstructs
    the MOFs into the various building units and creates a three json files

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
    ase_atoms_dic = {}
    search_data1 = {}
    search_data2 = {}
    seen = []
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)
        # porosity_dic = {}
    else:
        try:
            ase_atoms_dic = read_write.load_data(
                result_folder+'/ase_atoms_building_units.json')
            search_data1 = read_write.load_data(
                result_folder+'/sbus_and_linkers.json')
            search_data2 = read_write.load_data(
                result_folder+'/cluster_and_ligands.json')
            metal_info = read_write.load_data(result_folder+'/metal_info.json')
            seen = list(search_data1.keys())
        except Exception:
            pass

    for cif_file in cif_files:
        try:
            tmp_metal = {}
            encoder = read_write.AtomsEncoder
            base_name = cif_file[:cif_file.rindex('.')].split('/')[-1]
            ase_atom = read(cif_file)
            ase_atom = remove_guest(ase_atom)
            if not base_name in seen:
                data_to_json1, structural_data_1 = sbu_data(ase_atom)
                data_to_json2, structural_data_2 = ligand_data(ase_atom)
                search_data1[base_name] = data_to_json1
                search_data2[base_name] = data_to_json2
                structural_data = dict(merge_two_dicts(
                    structural_data_1, structural_data_2))
                ase_atoms_dic[base_name] = structural_data

                read_write.append_json_atom(ase_atoms_dic, encoder,
                                            result_folder+'/ase_atoms_building_units.json')
                read_write.append_json(
                    search_data1,  result_folder+'/sbus_and_linkers.json')
                read_write.append_json(
                    search_data2, result_folder+'/cluster_and_ligands.json')
        except Exception:
            pass

    if verbose:
        print(f"Saved results to {result_folder}")
    return


def main():
    '''
    CLI for deconstructing MOFs to building units.
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
