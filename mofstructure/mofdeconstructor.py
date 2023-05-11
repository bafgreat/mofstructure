
"""
Function for deconstructing MOFs into building units
There are three types of building units
1) The organic linkers, which contains the all atoms of the organic ligands
2) The metal sbu, which contains the metal cluster found in the MOF
3) The organic sbu, which is the fragment of the organic linker cut at the point of extension
"""
from ase.io import read
import numpy as np
from ase import neighborlist, geometry
from ase.data import chemical_symbols, covalent_radii, atomic_numbers
from scipy import sparse
from openbabel import pybel as pb
from openbabel import openbabel as ob
from rdkit.Chem import rdmolfiles, inchi, rdDetermineBonds
from rdkit import Chem
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.local_env import JmolNN
from math import sin, cos, pi

__name__ = "MOF.structure"
__author__ = "Dinga Wonanke"


def transitionMetals():
    '''
    Function containing list of symbols of all metals found in the periodic table
    '''
    metal = [symbol for symbol in chemical_symbols if symbol not in [chemical_symbols[main_index]
                                                                     for main_index in [1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 33, 34, 35, 36, 52, 53, 54, 85, 86]]]
    return metal

def covalent_radius(element):
    '''
    Extract ase covelent radii for an element
    '''
    a_n = atomic_numbers[element]
    return covalent_radii[a_n]

def ase2xyz(atoms):
    """
    Create and xyz string from an ase atom object to be compatible with
    pybel in order to perfom some cheminformatics.

    """

    if any(atoms.get_pbc()):
        raise RuntimeError(" Does not support periodic systems!")
    num_atoms = len(atoms)
    elements = atoms.get_chemical_symbols()
    all_atoms = zip(elements, atoms.get_positions())
    a_str = str(num_atoms) + "\n" + "\n"
    for atom in all_atoms:
        a_str += atom[0] + " " + " ".join([str(x) for x in atom[1]]) + "\n"
    return a_str[:-1]


def obmol2rdkit(OBmol):
    '''
    A simple function to convert openbabel mol to RDkit mol
    '''
    obconverted = ob.OBConversion()
    obconverted.SetOutFormat('sdf')
    outmdl = obconverted.WriteString(OBmol)
    rdmol = Chem.rdmolfiles.MolFromMolBlock(outmdl)
    return rdmol

def compute_inchis(obmol):
    '''
    Openbabel function to extract inchi and inchikey for molecules
    '''
    conv = ob.OBConversion()
    conv.SetOutFormat("inchi")
    inchi = conv.WriteString(obmol).rstrip()
    inchi = inchi.split('InChI=')[1]
    conv.SetOptions("K", conv.OUTOPTIONS)
    inchikey = conv.WriteString(obmol).rstrip()
    return inchi,  inchikey

def compute_smi(obmol):
    '''
    Openbabel function to extract smile strings for molecules
    '''
    conv = ob.OBConversion()
    conv.SetOutFormat("smi")
    smi = conv.WriteString(obmol).rstrip()
    return smi

def ase2pybel(atoms):
    """
    As simple script to convert from ase atom object to pybel
    Parameters
    ----------
    atoms : ASE atoms object

    Returns
    -------
    pybel: pybel molecular object.
    """
    a_str = ase2xyz(atoms)
    pybel_mol = pb.readstring("xyz", a_str)

    return pybel_mol

def max_index(lists):
    for i, elt in enumerate(lists):
        if elt == max(lists):
            return i

def compute_openbabel_cheminformatic(ase_atom):
    '''
    A procedure that uses openbabel to generate
    computer readable file formats
    '''

    new_ase_atom = wrap_systems_in_unit_cell(ase_atom)

    new_ase_atom.set_pbc(False)
    pybel_mol = ase2pybel(new_ase_atom)

    obmol = pybel_mol.OBMol
    inchi,  inchikey = compute_inchis(obmol)
    smi = compute_smi(obmol)
    return smi, inchi, inchikey

def compute_rdkit_cheminformatic(ase_atom):
    '''
    A function that converts and ase atom into rkdit mol to extract some cheminformatic 
    information.The script begins by taking an ase atom an creating an xyz string 
    of the molecule, which is stored to memory. The xyz string is then read using 
    the MolFromXYZBlock function in rdkit, which reads an xyz string of the molecule.
    However, this does not include any of the bonding information. To do this the following 
    commands are parsed to the rdkit mol,bond_moll = Chem.Mol(rdkit_mol)
    rdDetermineBonds.DetermineBonds(bond_moll,charge=0)
    https://github.com/rdkit/UGM_2022/blob/main/Notebooks/Landrum_WhatsNew.ipynb

    Parameters:
    -----------
    ASE atoms

    Returns
    -------
    smile string, inchi, inchikey
    '''
    if any(ase_atom.get_pbc()):
        ase_atom.set_pbc(False)

    xyz_string = ase2xyz(ase_atom)
    rdkit_xyz = Chem.MolFromXYZBlock(xyz_string)
    rdkit_mol = Chem.Mol(rdkit_xyz)
    rdDetermineBonds.DetermineConnectivity(rdkit_mol)
    rdDetermineBonds.DetermineBonds(rdkit_mol, charge=0)
    # Remove all hydrogens in the smile string
    smi = Chem.MolToSmiles(Chem.RemoveHs(rdkit_mol))
    InChi = inchi.MolToInchi(rdkit_mol)
    InChiKey = inchi.MolToInchiKey(rdkit_mol)
    return smi, InChi, InChiKey

def compute_ase_neighbour(ase_atom):
    '''
    Create a connectivity graph using ASE neigbour list.

    Parameters:
    -----------
    ASE atoms

    Returns
    -------
    Returns a python dictionary, wherein each atom index is key and the value
    are the indices of it neigbours.
    e.g.
    atom_neighbors ={0:[1,2,3,4], 1:[3,4,5]...}
    '''
    atom_neighbors = {}
    cutOff = neighborlist.natural_cutoffs(ase_atom)

    neighborList = neighborlist.NeighborList(cutOff,
                                             self_interaction=False,
                                             bothways=True)
    neighborList.update(ase_atom)
    matrix = neighborList.get_connectivity_matrix(sparse=False)

    for atoms in ase_atom:
        connectivity, _ = neighborList.get_neighbors(atoms.index)
        atom_neighbors[atoms.index] = connectivity

    return atom_neighbors, matrix

def matrix2dict(bond_matrix):
    '''
    A simple procedure to convert an adjacent matrix to a python dictionary
    Parameters:
    -----------
    Bond matrix
    type: nxn ndarray

    Returns
    -------
    python dictionary
    '''
    graph = {}
    for idx, row in enumerate(bond_matrix):
        temp = []
        for r in range(len(row)):
            if row[r] != 0:
                temp.append(r)
        graph[idx] = temp
    return graph

def dfsutil_graph_method(graph, temp, node, visited):
    '''
    Depth-first search graph algorithm for traversing graph data structures.
    I starts at the root 'node' and explores as far as possible along
    each branch before backtracking.
    It is used here a a util for searching connected components in the MOF graph
    Parameters:
    -----------
    graph: any python dictionary
    temp: a python list to hold nodes that have been visited
    node: a key in the python dictionary (graph), which is used as the starting or root node
    visited: python list containing nodes that have been traversed.
    Returns
    -------
    python dictionary
    '''
    visited[node] = True
    temp.append(node)
    for i in graph[node]:
        if visited[i] == False:
            temp = dfsutil_graph_method(graph, temp, i, visited)
    return temp

def remove_unbound_guest(ase_atom):
    '''
    A simple script to remove guest from a metal organic framework.
    1)It begins by computing a connected graph component of all the fragments in the system
    using ASE neighbour list.
    2) Secondly it selects indicies of connected components which contain a metal
    3)if the there are two or more components, we create a pytmagen graph for each components and filter out all components that are not polymeric
    4) If there are two or more polymeric components, we check wether these systems there are identical or different
    and select only unique polymeric components
    '''

    atom_neighbors, _ = compute_ase_neighbour(ase_atom)

    fragments = ConnectedComponents(atom_neighbors)

    if len(fragments) == 1:
        return [atom.index for atom in ase_atom]
    else:

        chemical_symbols = [[ase_atom[i].symbol for i in fragment if ase_atom[i].symbol in transitionMetals()]
                            for fragment in fragments]
        atomic_indices = [i for i in range(
            len(chemical_symbols)) if len(chemical_symbols[i]) > 0]

        polymeric_indices = []
        if len(atomic_indices) > 1:
            for i in atomic_indices:
                super_cell = ase_atom[fragments[i]]*(2, 1, 1)
                pymat_graph = ConnectedComponents(
                    compute_ase_neighbour(super_cell))
                if len(pymat_graph) == 0:
                    polymeric_indices.append(i)

            if len(polymeric_indices) > 1:
                Graphs = [StructureGraph.with_local_env_strategy(AseAtomsAdaptor.get_structure(
                    ase_atom[fragments[i]]), JmolNN()) for i in polymeric_indices]
                temp_indices = [polymeric_indices[0]]
                unique = [Graphs[0]]
                for k in range(len(polymeric_indices[1:])):
                    graph = Graphs[k]
                    if True in [graph.diff(j) for j in unique]:
                        unique.append(Graphs[k])
                        temp_indices.append(polymeric_indices[i])
                polymeric_indices = temp_indices
        else:
            polymeric_indices = atomic_indices
        mof_indices = []
        for i in polymeric_indices:
            mof_indices.extend(fragments[i])

        return mof_indices

def ConnectedComponents(graph):
    '''
    Find the connected fragments in a graph. Should work for any graph defined as a dictionary
    Parameters:
    -----------
    A graph in the form of dictionary
    e.g.
    graph = {1:[0,1,3], 2:[2,4,5]}

    Returns
    -------
    Returns a python list of list of connected components
    These correspond to individual molecular fragments.
    cc = [[1,2],[1,3,4]]
    '''
    visited = []
    cc = []
    for i in list(graph.keys()):
        visited.append(False)
    for v in list(graph.keys()):
        if visited[v] == False:
            temp = []
            cc.append(dfsutil_graph_method(graph, temp, v, visited))
    return cc

def check_planarity(p1, p2, p3, p4):
    '''
    A simple procedure to check whether a point is planar to three other points.
    Important to distinguish porphyrin type metals
    Parameters:
    -----------
    p1, p2, p3, p4 : ndarray containing x,y,z values

    Returns
    -------
    Boolean
    True: planar
    False: noneplanar
    '''
    planar = False
    a1 = p2[0] - p1[0]
    b1 = p2[1] - p1[1]
    c1 = p2[2] - p1[2]
    a2 = p3[0] - p1[0]
    b2 = p3[1] - p1[1]
    c2 = p3[2] - p1[2]
    # ------------------------------------------------
    a = b1*c2 - b2*c1
    b = a2*c1 - a1*c2
    c = a1*b2 - b1*a2
    d = round((-a*p1[0] - b*p1[1] - c*p1[2]), 0)
    # ------------------------------------------------
    factor = round((a*p4[0] + b*p4[1] + c*p4[2]), 0)
    verify = factor + d
    # ------------------------------------------------
    if verify == 0:
        planar = True
    return planar

def metal_in_porphyrin(ase_atom, graph):
    '''
    https://en.wikipedia.org/wiki/Transition_metal_porphyrin_complexes
    Check whether a metal is found at the centre of a phorphirin

    Parameters:
    -----------
    ase_atom: ASE atom
    graph: python dictionary containing neigbours

    Returns
    -------
    list of indices consisting of index of metal atoms found in the ASE atom
    '''
    all_porphyrin = []
    all_metal_symbols = [
        atom.index for atom in ase_atom if atom.symbol in transitionMetals()]
    for idx in all_metal_symbols:
        connected = graph[idx]
        all_nitrogens = [i for i in connected if ase_atom[i].symbol == 'N']
        if len(all_nitrogens) == 4:
            p1, p2, p3, p4 = ase_atom[all_nitrogens[0]].position, ase_atom[all_nitrogens[1]
                                                                 ].position, ase_atom[all_nitrogens[2]].position, ase_atom[all_nitrogens[3]].position
            planarity = check_planarity(p1, p2, p3, p4)
            if planarity:
                all_porphyrin.append(idx)
                all_porphyrin.extend(all_nitrogens)

    return all_porphyrin

def metal_in_porphyrin2(ase_atom, graph):
    '''
    https://en.wikipedia.org/wiki/Transition_metal_porphyrin_complexes
    Check whether a metal is found at the centre of a phorphirin

    Parameters:
    -----------
    ase_atom: ASE atom
    graph: python dictionary containing neigbours

    Returns
    -------
    list of indices consisting of index of metal atoms found in the ASE atom
    '''
    old_CC = ConnectedComponents(graph)
    all_porphyrin = []
    all_metal_symbols = [
        atom.index for atom in ase_atom if atom.symbol in transitionMetals()]
    metal_tmp = []
    N_tmp = []
    for idx in all_metal_symbols:
        connected = graph[idx]
        all_nitrogens = [i for i in connected if ase_atom[i].symbol == 'N']
        if len(all_nitrogens) == 4:
            metal_tmp.append(idx)
            N_tmp.extend(all_nitrogens)
    new_atom_indices = [i.index for i in ase_atomif i.index not in metal_tmp]
    tmp_atom = ase_atom[new_atom_indices]
    atom_neighbors, _ = compute_ase_neighbour(tmp_atom)
    CC = ConnectedComponents(atom_neighbors)
    if len(CC) == len(old_CC):
        all_porphyrin = metal_tmp + N_tmp
    return all_porphyrin

def move2front(index_value, coords):
    '''
    Parameters:
    -----------
    index_value: index of item to move to the from
    coords: list of coords

    Returns
    -------
    Move an index from any position in the list to the front
    The function is important to set the cell of a rodmof to point in the 
    a-axis. Such that the system can be grow along this axis
    '''
    if any(isinstance(el, list) for el in coords):
        for data in coords:
            data.insert(0, data.pop(index_value))
    else:
        coords.insert(0, coords.pop(index_value))
    return coords

def mof_regions(ase_atom, CC, Removed_dict):
    all_regions = {}
    all_pm_structures = [sorted(ase_atom[i].symbols) for i in CC]
    for i in range(len(all_pm_structures)):
        temp = []
        for j in range(len(all_pm_structures)):
            if all_pm_structures[i] == all_pm_structures[j]:
                temp.append(j)
        if not temp in all_regions .values():
            all_regions [i] = temp
    Xis_regions = {}
    for idx in range(len(all_regions .keys())):
        frag = list(all_regions .keys())[idx]
        components = CC[all_regions [frag][0]]
        Xis = [[Removed_dict[j] for j in Removed_dict.keys() if j in comp]
               for comp in components]
        Xis_regions[idx] = Xis
    return all_regions , Xis_regions

def find_carboxylates(ase_atom):
    '''
    A simple aglorimth to search for carboxylates found in the system.
    Parameters:
    -----------
    ase_atom: ASE atom

    Returns
    -------
    dictionary of key = carbon index and values = oxygen index
    '''
    graph, _ = compute_ase_neighbour(ase_atom)
    carboxyl = {}
    for atoms in ase_atom:
        if atoms.symbol == 'C':
            index = atoms.index
            oxygen = [i for i in graph[index] if ase_atom[i].symbol == 'O']
            if len(oxygen) == 2:
                oxy_metal = sum(
                    [[j for j in graph[i] if ase_atom[j].symbol in transitionMetals()] for i in oxygen], [])
                if len(oxy_metal) > 0:
                    carboxyl[index] = oxygen
    return carboxyl

def find_carbonyl_sulphate(ase_atom):
    '''
    A simple aglorimth to search for Carbonyl sulphate  found in the system.
       O
       |
    -C-S
       |
       O
    Parameters:
    -----------
    ase_atom: ASE atom

    Returns
    -------
    dictionary of key = carbon index and values = oxygen index
    '''
    graph, _ = compute_ase_neighbour(ase_atom)
    sulphate = {}
    for atoms in ase_atom:
        if atoms.symbol == 'S':
            index = atoms.index
            oxygen = [i for i in graph[index] if ase_atom[i].symbol == 'O']
            carbon = [i for i in graph[index] if ase_atom[i].symbol == 'C']
            if len(oxygen) >= 1:
                oxy_metal = sum(
                    [[j for j in graph[i] if ase_atom[j].symbol in transitionMetals()] for i in oxygen], [])
                if len(oxy_metal) > 0:
                    if len(carbon) > 0:
                        sulphate[index] = oxygen
    return sulphate

def secondary_building_units(ase_atom):
    """
    1) Search for all carboxylate that are connected to a metal.
       Cut at the position between the carboxylate carbon and the adjecent carbon.
    2) Find all Nitrogen connected to metal. Check whether the nitrogen is in the
       centre of a porphirin ring. If no, cut at nitrogen metal bond.
    3) Look for oxygen that is connected to metal and two carbon. cut at metal oxygen bond
    Parameters:
    -----------
    ase_atom: ASE atom

    Returns
    -------
     CC : list of connected components, in which each list contains atom indices
     Removed_dict : Dictionary containing point of disconnection
     Porphyrin_checker : Boolean showing whether the metal is in the centre of a porpherin
     Regions : Dictionary of regions.
    """
    graph, bond_matrix = compute_ase_neighbour(ase_atom)
    porphyrin_checker = metal_in_porphyrin2(ase_atom,  graph)
    Removed_dict = {}
    all_regions  = {}
    To_remove = []
    carboxylates = find_carboxylates(ase_atom)
    Sulphate = find_carbonyl_sulphate(ase_atom)
    for atoms in graph:
        if atoms in list(carboxylates.keys()):
            '''
            Remove carboxylate
            '''
            connected = graph[atoms]
            car_indx = [i for i in connected if ase_atom[i].symbol == 'C']
            N_indx = [i for i in connected if ase_atom[i].symbol == 'N']
            S_indx = [i for i in connected if ase_atom[i].symbol == 'S']
            if len(car_indx) == 1:
                To_remove.append([atoms] + car_indx)
                Removed_dict[atoms] = car_indx[0]
            if len(N_indx) == 1:
                To_remove.append([atoms] + N_indx)
                Removed_dict[atoms] = N_indx[0]
            if len(S_indx) == 1:
                To_remove.append([atoms] + S_indx)
                Removed_dict[atoms] = S_indx[0]

        if ase_atom[atoms].symbol == 'C':
            '''
            Remove carbon that is connected to one oxygen which is connected to a metal
                  R
                  |
                R-C=O
             '''
            connected = graph[atoms]
            oxygens = [i for i in connected if ase_atom[i].symbol == 'O']
            if len(oxygens) == 1:
                oxy_metal = [
                    i for i in oxygens if ase_atom[i].symbol in transitionMetals()]
                if len(oxy_metal) == 1:
                    Removed_dict[oxygens[0]] = oxy_metal[0]
                    To_remove.append([oxygens] + oxy_metal)

        if atoms in list(Sulphate.keys()):
            '''
            Remove sulphate
            '''
            connected = graph[atoms]
            car_indx = [i for i in connected if ase_atom[i].symbol == 'C']
            oxygen = Sulphate[atoms]
            if len(car_indx) == 1:
                Removed_dict[atoms] = car_indx[0]
                To_remove.append([atoms] + car_indx)
            if len(car_indx) > 1:
                for oxy in oxygen:
                    metal = [i for i in graph[oxy]
                             if ase_atom[i].symbol in transitionMetals()]
                    if len(metal) > 0:
                        for met in metal:
                            Removed_dict[oxy] = met
                            To_remove.append([oxy] + [met])

        if ase_atom[atoms].symbol == 'O':
            seen = sum(list(carboxylates.values())+list(Sulphate.values()), [])
            if not atoms in seen:
                connected = graph[atoms]
                metal = [
                    i for i in connected if ase_atom[i].symbol in transitionMetals()]
                Nitrogen = [i for i in connected if ase_atom[i].symbol == 'N']
                carbon = [i for i in connected if ase_atom[i].symbol ==
                          'C' and i not in list(carboxylates.keys())]
                if len(metal) == 1 and len(carbon) == 1:
                    Removed_dict[atoms] = carbon[0]
                    To_remove.append([atoms] + carbon)
                if len(metal) == 1 and len(Nitrogen) == 1:
                    n_carbon = [i for i in graph[Nitrogen[0]]
                                if ase_atom[i].symbol == 'C' and i not in list(carboxylates.keys())]
                    n_nitrogen = [i for i in graph[Nitrogen[0]]
                                  if ase_atom[i].symbol == 'N']
                    n_sulphur = [i for i in graph[Nitrogen[0]]
                                 if ase_atom[i].symbol == 'S' and i not in list(Sulphate.keys())]
                    if len(n_carbon) > 1:
                        Removed_dict[atoms] = Nitrogen[0]
                        To_remove.append([atoms] + Nitrogen)

                    elif len(n_nitrogen) > 1:
                        Removed_dict[atoms] = Nitrogen[0]
                        To_remove.append([atoms] + Nitrogen)
                    elif len(n_sulphur) > 1:
                        Removed_dict[atoms] = Nitrogen[0]
                        To_remove.append([atoms] + Nitrogen)

        if ase_atom[atoms].symbol == 'N':
            connected = graph[atoms]
            metal = [i for i in connected if ase_atom[i].symbol in transitionMetals()]
            if atoms not in porphyrin_checker:
                if len(metal) == 1:
                    Removed_dict[atoms] = metal[0]
                    To_remove.append([atoms] + metal)

        if ase_atom[atoms].symbol == 'S':
            connected = graph[atoms]
            metal = [i for i in connected if ase_atom[i].symbol in transitionMetals()]
            if len(metal) > 0:
                for met in metal:
                    Removed_dict[atoms] = met
                    To_remove.append([atoms, met])

        if ase_atom[atoms].symbol == 'P':
            '''
            Find the carbon closest to P, which is not bonded to a metal and cut
            1) Look for phosphorous atoms
            2) Look for all it's neigbours
            3) Look for neigbours that are not connected to metal or hydrogen.
            '''
            connected = graph[atoms]
            print(connected)
            not_connected_to_metal_or_hygrogen = [
                [i for i in graph[j] if ase_atom[i].symbol not in transitionMetals() or ase_atom[i].symbol != 'H'] for j in connected]

            metal_oxy = [[i for i in graph[j] if ase_atom[i].symbol in transitionMetals()]
                         for j in connected]

            metal = sum(metal_oxy, [])
            closest_atoms = sum(
                [[i for i in graph[j] if i != atoms and not ase_atom[i].symbol in transitionMetals()] for j in connected], [])

            if len(metal) > 0:
                car_indx = sum([[i for i in graph[j] if i in connected]
                               for j in closest_atoms], [])
                for frag in car_indx:
                    Removed_dict[atoms] = frag
                    Removed_dict[frag] = atoms
                    To_remove.append([atoms, frag])

        if ase_atom[atoms].symbol == 'B':
            '''
            Find the carbon closest to P, which is  bonded to a metal and cut
            '''
            connected = [i for i in graph[atoms]
                         if ase_atom[i].symbol not in transitionMetals()]
            metal_oxy = [[[i, j] for i in graph[j] if ase_atom[i].symbol in transitionMetals()]
                         for j in connected]

            metal_connect = sum(metal_oxy, [])
            if len(metal_connect) > 0:
                for frag in metal_connect:
                    Removed_dict[frag[0]] = frag[1]
                    To_remove.append(frag)

    for bonds in To_remove:
        bond_matrix[bonds[0], bonds[1]] = 0
        bond_matrix[bonds[1], bonds[0]] = 0

    new_ase_graph= matrix2dict(bond_matrix)
    try:
        CC = ConnectedComponents(new_ase_graph)
    except:
        import networkx as nx
        N_Graph = nx.from_dict_of_lists(new_ase_graph)
        CC = [list(i) for i in list(nx.connected_components(N_Graph))]

    all_pm_structures = [sorted(ase_atom[i].symbols) for i in CC]
    for i in range(len(all_pm_structures)):
        temp = []
        for j in range(len(all_pm_structures)):
            if all_pm_structures[i] == all_pm_structures[j]:
                temp.append(j)
        if not temp in all_regions .values():
            all_regions [i] = temp

    return CC, Removed_dict, porphyrin_checker , all_regions 

def ligands_and_metal_clusters(ase_atom):
    '''
    Start by checking whether there are more than 2 layers
    if yes, select one
    Here we select the largest connected component
    '''
    graph, bond_matrix = compute_ase_neighbour(ase_atom)

    porphyrin_checker = metal_in_porphyrin2(ase_atom,  graph)

    all_regions  = {}
    Removed_dict = {}
    To_remove = []
    carboxylates = find_carboxylates(ase_atom)
    Sulphate = find_carbonyl_sulphate(ase_atom)
    Seen_oxygen = []
    for atoms in graph:
        if atoms in list(carboxylates.keys()):
            '''
            Remove carboxylate
            '''
            oxygen = carboxylates[atoms]
            for oxy in oxygen:

                metal = [i for i in graph[oxy]
                         if ase_atom[i].symbol in transitionMetals()]
                if len(metal) > 0:
                    Seen_oxygen.append(oxy)
                    for met in metal:
                        To_remove.append([atoms] + [met])
                        Removed_dict[atoms] = met

        if atoms in list(Sulphate.keys()):
            '''
            Remove carboxylate
            '''
            oxygen = Sulphate[atoms]
            connected = graph[atoms]
            car_indx = [i for i in connected if ase_atom[i].symbol == 'C']
            for oxy in oxygen:
                metal = [i for i in graph[oxy]
                         if ase_atom[i].symbol in transitionMetals()]
                if len(metal) > 0:
                    for met in metal:
                        To_remove.append([oxy] + [met])
                        Removed_dict[oxy] = met

        if ase_atom[atoms].symbol == 'N':
            connected = graph[atoms]
            metal = [i for i in connected if ase_atom[i].symbol in transitionMetals()]
            if len(metal) > 0 and atoms not in porphyrin_checker:
                for met in metal:
                    To_remove.append([atoms, met])
                    Removed_dict[atoms] = met
                # Removed_dict[metal[0]] = atoms

        if ase_atom[atoms].symbol == 'S':
            connected = graph[atoms]
            metal = [i for i in connected if ase_atom[i].symbol in transitionMetals()]
            if len(metal) > 0:
                if len(metal) > 0:
                    for met in metal:
                        To_remove.append([atoms] + [met])
                        Removed_dict[atoms] = met

        if ase_atom[atoms].symbol == 'O':
            if not atoms in Seen_oxygen:
                connected = graph[atoms]
                metal = [
                    i for i in connected if ase_atom[i].symbol in transitionMetals()]
                Nitrogen = [i for i in connected if ase_atom[i].symbol == 'N']
                carbon = [i for i in connected if ase_atom[i].symbol == 'C']
                if len(metal) > 0 and len(carbon) == 1:
                    for met in metal:
                        Removed_dict[atoms] = met
                        To_remove.append([atoms, met])
                if len(metal) > 0 and len(Nitrogen) == 1:
                    n_carbon = [i for i in graph[Nitrogen[0]]
                                if ase_atom[i].symbol in ['C', 'S', 'N']]
                    if len(n_carbon) > 1:
                        for met in metal:
                            Removed_dict[atoms] = met
                        To_remove.append([atoms, met])

                if len(carbon) > 0:
                    if len(metal) > 0:
                        for met in metal:
                            To_remove.append([atoms] + [met])
                            Removed_dict[atoms] = met

        if ase_atom[atoms].symbol == 'P':
            '''
            Find the carbon closest to P, which is not bonded to a metal and cut
            '''
            connected = [i for i in graph[atoms]
                         if ase_atom[i].symbol not in transitionMetals()]
            metal_oxy = [[i for i in graph[j] if ase_atom[i].symbol in transitionMetals()]
                         for j in connected]

            metal = sum(metal_oxy, [])
            closest_atoms = sum([[[i, j] for i in graph[j] if i !=
                                atoms and ase_atom[i].symbol in transitionMetals()] for j in connected], [])
            for frag in closest_atoms:
                To_remove.append(frag)
                Removed_dict[frag[0]] = frag[1]

    for bonds in To_remove:
        bond_matrix[bonds[0], bonds[1]] = 0
        bond_matrix[bonds[1], bonds[0]] = 0

    new_ase_graph= matrix2dict(bond_matrix)
    try:
        CC = ConnectedComponents(new_ase_graph)
    except:
        import networkx as nx
        N_Graph = nx.from_dict_of_lists(new_ase_graph)
        CC = [list(i) for i in list(nx.connected_components(N_Graph))]

    all_pm_structures = [sorted(ase_atom[i].symbols) for i in CC]
    for i in range(len(all_pm_structures)):
        temp = []
        for j in range(len(all_pm_structures)):
            if all_pm_structures[i] == all_pm_structures[j]:
                temp.append(j)
        if not temp in all_regions .values():
            all_regions [i] = temp

    return CC, Removed_dict, porphyrin_checker, all_regions 


def Is_rod(metal_sbu):
    '''
    Simple test to check whether a metal sbu is a rodlike MOF
    '''
    Rod_check = []
    cells = [(2, 1, 1), (1, 2, 1), (1, 1, 2)]
    graph, _ = compute_ase_neighbour(metal_sbu)
    for index, ijk in enumerate(cells):
        rod = metal_sbu*ijk
        graph, _ = compute_ase_neighbour(rod)
        CC = ConnectedComponents(graph)
        if len(CC) == 1:
            Rod_check.append(index)
    return Rod_check

def is_ferrocene(metal_sbu, graph):
    '''
    A simple script to check whether a metal_sbu is ferrocene
    '''
    Check = []
    verdict = False
    All_connectivity = list(graph.values())
    for connectivity in All_connectivity:
        if len(connectivity) >= 10:
            carbons = 0
            for bondedAtomIndex in connectivity:
                if metal_sbu[bondedAtomIndex].symbol == 'C':
                    carbons += 1
            if carbons >= 10:
                verdict = True
            Check.append(verdict)
        Correct = False
        if True in Check:
            Correct = True
    return Correct


def is_paddlewheel(metal_sbu, graph):
    """
    Returns True if the atom is part of a paddlewheel motif
    """
    Check = []
    verdict = False
    All_connectivity = list(graph.values())
    for connectivity in All_connectivity:
        metalNeighbours = 0
        oxygenNeighbours = 0
        for bondedAtomIndex in connectivity:
            if metal_sbu[bondedAtomIndex].symbol in transitionMetals():
                metalNeighbours += 1
            if metal_sbu[bondedAtomIndex].symbol == 'O':
                oxygenNeighbours += 1
        if metalNeighbours == 1 and oxygenNeighbours == 4 and len(connectivity) >= 5 and len(connectivity) <= 6:
            verdict = True
        Check.append(verdict)
    Correct = False
    if True in Check:
        Correct = True
    return Correct

def is_Paddle_water(ase_atom, graph):
    """
    Returns True if the atom is part of a paddle wheel with water motif
    """
    Check = []
    metal = []
    verdict = False
    for atoms in ase_atom:
        if atoms.symbol in transitionMetals():
            index = atoms.index
            metal.append(index)
            connectivity = graph[index]
            if len([ase_atom[i].symbol for i in connectivity if ase_atom[i].symbol == 'O']) == 5:
                verdict = True
                Check.append(verdict)
    Correct = False
    if True in Check and len(metal) == 2:
        Correct = True
    return Correct

def is_UIO66(ase_atom, graph):
    """
    Returns True if the atom is part of a UIO66 motif
    """
    Check = []
    verdict = False
    All_connectivity = list(graph.values())
    for connectivity in All_connectivity:
        metalNeighbours = 0
        oxygenNeighbours = 0
        for bondedAtomIndex in connectivity:
            if ase_atom[bondedAtomIndex].symbol in transitionMetals():
                metalNeighbours += 1
            if ase_atom[bondedAtomIndex].symbol == 'O':
                oxygenNeighbours += 1
        if metalNeighbours == 4 and (oxygenNeighbours == 6 or oxygenNeighbours == 8) and (len(connectivity) == 10 or len(connectivity) == 12):
            verdict = True
        Check.append(verdict)
        Correct = False
        if True in Check:
            Correct = True
    return Correct

def is_IRMOF(ase_atom, graph):
    """
    Returns True if the atom is part of a IRMOF motif
    """
    Check = []
    verdict = False
    for atoms in ase_atom:
        if atoms.symbol == 'O':
            index = atoms.index
            connectivity = graph[index]
            if len(connectivity) == 4 and len([ase_atom[i].symbol for i in connectivity if ase_atom[i].symbol in transitionMetals()]) == 4:

                verdict = True
                Check.append(verdict)
    Correct = False
    if True in Check:
        Correct = True
    return Correct

def is_MOF32(ase_atom, graph):
    """
    Returns True if the atom is part of a MOF32 motif
    """
    Check = []
    metal = []
    verdict = False
    for atoms in ase_atom:
        if atoms.symbol in transitionMetals():
            index = atoms.index
            metal.append(index)
            connectivity = graph[index]
            if len([ase_atom[i].symbol for i in connectivity if ase_atom[i].symbol == 'O']) == 8:
                verdict = True
                Check.append(verdict)
    Correct = False
    if True in Check and len(metal) == 1:
        Correct = True
    return Correct

def Rod_Manipulation(ase_atom, checker):
    '''
    Script to adjust Rodlike sbus.
    1) Its collects the axis responsible for expanding the rod
    2) It shifts all coordinates to the axis
    3) It rotates the rod to lie in the directions of expansion,
    '''

    cell = ase_atom.get_cell().tolist()
    value_index = checker[0]
    cell_value = cell[value_index]
    new_cell = move2front(checker, cell_value)
    coords = ase_atom.positions - cell_value
    new_position = move2front(checker, coords.tolist())
    ase_atom.positions = new_position
    return ase_atom, new_cell

def find_unique_building_units(CC, Removed_dict, ase_atom, porphyrin_checker , all_regions,wrap_system = True ):
    '''
    Find Unique components
    Returns a list of unique molecules
    '''
    ASE_metal, ASE_linker = [], []
    building_unit_regions = {}
    for idx, key in enumerate(all_regions .keys()):
        frag = list(all_regions .keys())[idx]
        components = CC[all_regions [frag][0]]
        Molecule_to_write = ase_atom[components]
        if wrap_system:
            Molecule_to_write = wrap_systems_in_unit_cell(Molecule_to_write)
        smi, inchi, inchikey = compute_openbabel_cheminformatic(
            Molecule_to_write)
        Molecule_to_write.info['smi'] = smi
        Molecule_to_write.info['inchi'] = str(inchi)
        Molecule_to_write.info['inchikey'] = str(inchikey)
        Molecule_to_write.info['atom_indices_mapping'] = [
            CC[i] for i in all_regions [key]]
        metal = [i.index for i in Molecule_to_write if i.symbol in transitionMetals(
        ) and i.index not in porphyrin_checker]

        if len(metal) > 0:
            graph_sbu, _ = compute_ase_neighbour(Molecule_to_write)
            '''
            Check whether the metal sbu is a rod mof. If is it is rod mof,
            we rotate and aligne the sbu such that the axis of rotation will be the a-axis.
            '''
            if len(Is_rod(Molecule_to_write)) == 1:
                Molecule_to_write.info['sbu_type'] = 'is_rod'
       
            elif is_ferrocene(Molecule_to_write, graph_sbu):
                Molecule_to_write.info['sbu_type'] = 'is_ferrocene'
            elif is_paddlewheel(Molecule_to_write, graph_sbu):
                Molecule_to_write.info['sbu_type'] = 'is_paddlewheel'
            elif is_Paddle_water(Molecule_to_write, graph_sbu):
                Molecule_to_write.info['sbu_type'] = 'is_paddlewheel_with_water'
            elif is_UIO66(Molecule_to_write, graph_sbu):
                Molecule_to_write.info['sbu_type'] = 'is_UIO66'
            elif is_MOF32(Molecule_to_write, graph_sbu):
                Molecule_to_write.info['sbu_type'] = 'is_MOF32'
            elif is_IRMOF(Molecule_to_write, graph_sbu):
                Molecule_to_write.info['sbu_type'] = 'is_IRMOF'
            else:
                Molecule_to_write.info['sbu_type'] = 'still checking!'
            ASE_metal.append(Molecule_to_write)

        else:
            ASE_linker.append(Molecule_to_write)
        building_unit_regions[idx] = Molecule_to_write

    return ASE_metal, ASE_linker, building_unit_regions

def coordination_number(ase_atom):
    '''
    Extract coordination number of central metal
    '''
    CN = {}
    graph, _ = compute_ase_neighbour(ase_atom)
    metal_indices = [
        i.index for i in ase_atom if i.symbol in transitionMetals()]
    metal_elt = []
    CN = {}
    for i in metal_indices:
        if ase_atom[i].symbol not in metal_elt:
            metal_elt.append(ase_atom[i].symbol)
            CN[ase_atom[i].symbol] = len(graph[i])

    return metal_elt, CN

def wrap_systems_in_unit_cell(ase_atom, max_iter):
    '''
    A simple aglorithm to reconnnect all atoms wrapped in a periodic boundary condition such that all atoms outside the box will appear reconnected.
    '''
    new_position = geometry.wrap_positions(ase_atom.positions, ase_atom.cell, pbc=True, center=(
        0, 0, 0), pretty_translation=True, eps=1e-07)
    ase_atom.positions = new_position

    graph, bond_matrix = compute_ase_neighbour(ase_atom)

    for atom in graph:
        connected = graph[atom]
        for nl in connected:
            check = covalent_radius(
                ase_atom[atom].symbol) + covalent_radius(ase_atom[nl].symbol)+.3
            bond = round(ase_atom.get_distance(atom, nl), 2)
            if bond > check:
                bond_matrix[atom][nl] = 0

    new_ase_graph = matrix2dict(bond_matrix)
    CC = ConnectedComponents(new_ase_graph)
    iter = 0
    while len(CC) != 1:
        all_len = [len(i) for i in CC]
        max_index = all_len.index(max(all_len))
        Root = CC[max_index]
        CC.pop(max_index)
        All_sum = sum(CC, [])
        for atom in Root:
            connected = graph[atom]
            for nl in All_sum:
                if nl in connected:
                    v = ase_atom[nl].position - ase_atom[atom].position
                    mic_vector = geometry.find_mic(v, ase_atom.get_cell(), pbc=True)
                    ase_atom[nl].position = mic_vector[0] + ase_atom[atom].position

        graph, bond_matrix = compute_ase_neighbour(ase_atom)
        for atom in graph:
            connected = graph[atom]
            for nl in connected:
                check = covalent_radius(
                    ase_atom[atom].symbol) + covalent_radius(ase_atom[nl].symbol)+.3
                bond = round(ase_atom.get_distance(atom, nl), 2)
                if bond > check:
                    bond_matrix[atom][nl] = 0
        new_ase_graph = matrix2dict(bond_matrix)
        CC = ConnectedComponents(new_ase_graph)
        iter += 1
        if iter == max_iter:
            break

    return ase_atom

def remove_pbc(ase_atom):
    from scipy import sparse
    graph, bond_matrix = compute_ase_neighbour(ase_atom)
    n_components, component_list = sparse.csgraph.connected_components(
        bond_matrix)
    return
