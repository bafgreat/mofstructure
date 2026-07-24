#!/usr/bin/python
'''
Curation of covalent organic framework structures.

Wraps the deconstruction machinery for COF specific workflows, where the
building units are organic on both sides of the bond rather than metal and
linker, and collects the results into the tabular form used for building
datasets.
'''
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
import os
from ase.io import read
from functools import reduce
import operator
import argparse
import pandas as pd
import shutil
import tempfile
import json
import numpy as np
import logging
from typing import List, Tuple, Dict, Any
from ase.io import read
import mofstructure.filetyper as read_write
from mofstructure import mofdeconstructor


def find_CN_double_bonds(
    ase_atom,
    graph: np.ndarray,
    ) -> List[List[int]]:
    """
    A function to find Nitrogen connected by double bonds in carbon atoms in a given ASE Atoms object

            R
            |
        R-N=C-R

    **parameter:**
        ase_atom : ASE Atoms object
    **returns**
        bond_pairs of carbon and nitrogen atoms connected by double bonds
    """
    symbols = ase_atom.get_chemical_symbols()
    bond_pairs: List[List[int]] = []

    for c_idx, sym in enumerate(symbols):
        if sym != "C":
            continue

        neighbors = graph[c_idx]
        if len(neighbors) != 3:
            continue

        n_idx = -1
        n_count = 0
        h_count = 0

        for j in neighbors:
            sj = symbols[j]
            if sj == "N":
                n_count += 1
                n_idx = j
                if n_count > 1:  # early exit
                    break
            elif sj == "H":
                h_count += 1
                if h_count > 1:  # early exit
                    break

        if n_count == 1 and h_count == 1:
            bond_pairs.append([c_idx, n_idx])

    return bond_pairs

def secondary_building_units(ase_atom):
    """

    **parameters:**
        ase_atom: ASE atom

    **returns:**
        list_of_connected_components: list of connected components,in which each list contains atom indices
        bonds_to_break: List of lists of atom indices at breaking points
        porphyrin_checker: Boolean showing whether the metal is in the centre of a porpherin
        Regions: Dictionary of regions.
    """
    graph, bond_matrix = mofdeconstructor.get_neighbour_bond_matrix(ase_atom, skin=0.30, bo_step=0.10, aromatic=True)
    porphyrin_checker = mofdeconstructor.metal_in_porphyrin2(ase_atom, graph)
    bonds_to_break= []
    all_regions = {}
    cn_double_bonds = find_CN_double_bonds(ase_atom, graph)
    if len(cn_double_bonds) > 0:
        bonds_to_break.extend(cn_double_bonds)

    for bonds in bonds_to_break:
        bond_matrix[bonds[0], bonds[1]] = 0
        bond_matrix[bonds[1], bonds[0]] = 0

    new_ase_graph = mofdeconstructor.matrix2dict(bond_matrix)
    try:
        list_of_connected_components = mofdeconstructor.connected_components(new_ase_graph)
    except Exception:
        import networkx as nx
        N_Graph = nx.from_dict_of_lists(new_ase_graph)
        list_of_connected_components = [
            list(i) for i in list(nx.connected_components(N_Graph))]

    all_pm_structures = [sorted(ase_atom[i].symbols)
                         for i in list_of_connected_components]
    len_symbols = [len(i) for i in all_pm_structures]
    for i in range(len(all_pm_structures)):
        temp = []
        for j in range(len(all_pm_structures)):
            if all_pm_structures[i] == all_pm_structures[j]:
                temp.append(j)
        if temp not in all_regions .values():
            all_regions[i] = temp
    return [
        list_of_connected_components,
        bonds_to_break,
        porphyrin_checker, all_regions
    ]


def unique_building_units(list_of_connected_components,
                          bonds_to_break,
                          ase_atom,
                          porphyrin_checker,
                          all_regions,
                          wrap_system=True,
                          cheminfo=True,
                          add_dummy=False
                          ):
  _, cof_linkers, _ = mofdeconstructor.find_unique_building_units(list_of_connected_components,
                                                bonds_to_break,
                                                ase_atom,
                                                porphyrin_checker,
                                                all_regions,
                                                wrap_system=wrap_system,
                                                cheminfo=cheminfo,
                                                add_dummy=add_dummy)
  return cof_linkers




# ase_atom = read('CC1_P3_AA_kgm-refine-mono.cif')
# # graph, _ = mofdeconstructor.compute_ase_neighbour(ase_atom)

# # bonded_pair = find_CN_double_bonds(ase_atom, bonds_order)
# # print(atom_neighbors)
# # print(bonded_pair)

# cc, bonds_to_break, porphyrin_checker, all_regions = secondary_building_units(ase_atom)
# cof_linkers = unique_building_units(cc, bonds_to_break, ase_atom, porphyrin_checker, all_regions)
# for i, linker in enumerate(cof_linkers):
#     linker.write(f'test_{i}.xyz')

# print(all_regions)