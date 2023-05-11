
#!/usr/bin/python
"""
Create databse of building units
"""
import ase
from ase import Atoms
import json
import filetyper as File_typer
import glob
from ase.io import read
import mofdeconstructor as MOF_deconstructor
from ase.io import read, write
import sys
import os
__name__ = "MOF.structure"
__author__ = "Dinga Wonanke"


def put_contents(filename, output):
    with open(filename, 'w') as f:
        f.writelines(output)
    return


class AtomsEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, Atoms):
            coded = dict(positions=[list(pos) for pos in obj.get_positions()], cell=[
                         list(c) for c in obj.get_cell()], symbols=list(obj.get_chemical_symbols()))
            # coded['n_atoms'] = len(list(obj.get_chemical_symbols()))
            # coded['numbers'] = obj.get_atomic_numbers().tolist()

            keys = list(obj.info.keys())
            if 'atom_indices_mapping' in keys:
                info = obj.info
                coded.update(info)
            return coded
        if isinstance(obj, ase.spacegroup.Spacegroup):
            return obj.todict()
        return json.JSONEncoder.default(self, obj)


def sbu_data(experiment_atom):

    CC, Removed_dict, Porpyrin_checker, Regions = MOF_deconstructor.secondary_building_units(
        experiment_atom)

    ASE_metal, ASE_linker, building_unit_regions = MOF_deconstructor.find_unique_building_units(
        CC, Removed_dict, experiment_atom, Porpyrin_checker, Regions)
    for i in range(len(ASE_metal)):
        ASE_metal[i].write('metal_sbu_'+str(i+1)+'.xyz')
    for j in range(len(ASE_linker)):
        ASE_linker[j].write('organic_sbu_'+str(j+1)+'.xyz')

    return


def ligand_data(experiment_atom):

    CC, Removed_dict, Porpyrin_checker, Regions = MOF_deconstructor.ligands_and_metal_clusters(
        experiment_atom)

    ASE_metal, ASE_linker, building_unit_regions = MOF_deconstructor.find_unique_building_units(
        CC, Removed_dict, experiment_atom, Porpyrin_checker, Regions)
    for i in range(len(ASE_metal)):
        ASE_metal[i].write('metal_cluster_'+str(i+1)+'.xyz')
    for j in range(len(ASE_linker)):
        ASE_linker[j].write('organic_ligand_'+str(j+1)+'.xyz')
    return


experiment_atom = read('../data/EDUSIF.cif')
index_non_guest = MOF_deconstructor.remove_unbound_guest(experiment_atom)

experiment_atom[index_non_guest].write('../data/Test.cif')

# experiment_atom = read('Test.cif')
# sbu_data(experiment_atom)
# ligand_data(experiment_atom)
