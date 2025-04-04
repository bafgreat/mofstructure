#!/usr/bin/python
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
import numpy as np
from pymatgen.analysis.local_env import JmolNN
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.io.ase import AseAtomsAdaptor
from ase.data import chemical_symbols, covalent_radii, atomic_numbers
from ase import neighborlist, geometry
from ase.neighborlist import NeighborList
# from ase import Atoms

try:
    from openbabel import pybel as pb
    from openbabel import openbabel as ob
except ModuleNotFoundError:
    print('install openbabel if you wish to use compute_openbabel_cheminformatic function')
    print("pip install openbabel-wheel==3.1.1.16")
try:
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds
except ModuleNotFoundError:
    pass


def transition_metals():
    '''
    A function that returns the list of all chemical symbols of
    metallic elements present in the periodic table
    '''
    metal = [symbol for symbol in chemical_symbols if symbol not in [chemical_symbols[main_index]
                                                                     for main_index in [1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 33, 34, 35, 36, 52, 53, 54, 85, 86]]]
    return metal


def inter_atomic_distance_check(ase_atom):
    """
    A function that checks whether any pair of atoms (except those with exactly one hydrogen)
    has a distance below 0.90 Å. Only unique pairs are checked to avoid redundancy.

    **Parameters**
    ase_atom : ase.Atoms
        An ASE Atoms object containing atomic positions and chemical symbols.

    **Returns**
    bool
        Returns False if any applicable atom pair has a distance below 0.90 Å.
        Otherwise, returns True.
    """
    graph, _ = compute_ase_neighbour(ase_atom)
    for i, neighbours in graph.items():
        distances = ase_atom.get_distances(i, neighbours, mic=True)
        for j, d in zip(neighbours, distances):
            if i >= j:
                continue

            i_symbol = ase_atom[i].symbol
            j_symbol = ase_atom[j].symbol

            if (i_symbol == 'H') ^ (j_symbol == 'H'):
                continue

            if d < 0.9:
                return False
    return True


def inter_atomic_distance_check2(ase_atom):
    '''
    Checks whether any non-hydrogen atom in the ASE atoms object has a neighbor
    within 0.90 Å (ignoring self-distances). This is used to detect overly close contacts,
    while ignoring R–H bonds (i.e., cases where the central atom is hydrogen).

    **parameters:**
        ase_atom : ASE atoms object

    **returns:**
        bool: False if any non-hydrogen atom has a neighbor (other than itself) within 0.90 Å; True otherwise.
    '''
    distances = ase_atom.get_all_distances(mic=True)
    symbols = np.array(ase_atom.get_chemical_symbols())
    print(len(distances))

    non_h_mask = symbols != 'H'
    off_diag_mask = ~np.eye(len(symbols), dtype=bool)
    check_mask = non_h_mask[:, None] & off_diag_mask

    if np.any(distances[check_mask] < 0.90):
        return False
    return True

def inter_atomic_distance_iterative(ase_atom):
    '''
    A function that checks whether two atoms are within a distance 1.0 Amstrong unless it is an R-H bond

    **parameters:**
        ase_atom : ASE atoms object

    **returns**
        boolean :
    '''
    valid = True
    distances = ase_atom.get_all_distances(mic=True)
    for i in range(len(distances)):
        if ase_atom[i].symbol != 'H':
            for j in range(len(distances[i])):
                if i != j:
                    if distances[i, j] < 0.90:
                        valid = False
                        break
    return valid


def covalent_radius(element):
    '''
    A function that returns the covalent radius from the chemical symbol
    of an element.

    **parameters:**
        element: chemical symbol of an atom
        type.string

    **returns:**
        covalent_radii: covalent radius of the elements.
        type.float
    '''
    a_n = atomic_numbers[element]
    return covalent_radii[a_n]


def ase_2_xyz(atoms):
    """
    Create an xyz string from an ase atom object to be compatible with
    pybel in order to perfom some cheminformatics.

    **parameters:**
        atoms : ASE atoms object

    **returns:**
        a_str : string block of the atom object in xyz format

    """
    if any(atoms.get_pbc()):
        raise RuntimeError(" Does not support periodic systems!")
    num_of_atoms = len(atoms)
    all_chemical_symbols = atoms.get_chemical_symbols()
    all_atoms = zip(all_chemical_symbols, atoms.get_positions())
    a_str = str(num_of_atoms) + "\n" + "\n"
    for atom in all_atoms:
        a_str += atom[0] + " " + " ".join([str(x) for x in atom[1]]) + "\n"
    return a_str[:-1]


def obmol_2_rdkit(obmol):
    '''
    A simple function to convert openbabel mol to RDkit mol
    A function that takes an openbabel molecule object and
    converts it to an RDkit molecule object. The importance of this
    function lies in the fact that there is no direct were to convert from
    an ase atom type to an rdkit molecule object. Consequently, the easiest
    approach it is firt convert the system to an openbabel molecule object
    using the function ase_2_pybel to convert to the pybel molecule object.
    Then using the following code obmol = pybel_mol.OBMol to obtain the
    openbabel molecule object that can then be used to convert to the
    rdkit molecule object.

    **parameters:**
        obmol : openbabel molecule object

    **returns:**
        rdmol : rdkit molecule object.

    '''
    if Chem is None or rdDetermineBonds is None:
        print('install rdkit if you wish to use compute_cheminformatic_from_rdkit function')
    else:
        obconverted = ob.OBConversion()
        obconverted.SetOutFormat('sdf')
        outmdl = obconverted.WriteString(obmol)
        rdmol = Chem.MolFromMolBlock(outmdl)
        return rdmol

    return None


def compute_inchis(obmol):
    '''
    A function to compute the cheminformatic inchi and inchikey using
    openbabel. The inchi and inchikey are IUPAC cheminformatic identifiers
    of molecules. They are important to easily search for molecules in databases.
    For MOFs, this is quite important to search for different secondary building units
    (sbu) and ligands found in the MOF.
    More about inchi can be found in the following link

    https://iupac.org/who-we-are/divisions/division-details/inchi

    **parameters:**
        obmol: openbabel molecule object

    **returns:**
        inChi, inChiKey: iupac inchi and inchikeys for thesame molecule.
    '''
    conv = ob.OBConversion()
    conv.SetOutFormat("inchi")
    inChi = conv.WriteString(obmol).rstrip()
    inChi = inChi.split('InChI=')[1]
    conv.SetOptions("K", conv.OUTOPTIONS)
    inChiKey = conv.WriteString(obmol).rstrip()
    return inChi,  inChiKey


def compute_smi(obmol):
    '''
    A function to compute the SMILES (Simplified molecular-input line-entry system) notation
    of a molecule using openbabel. The SMILES is a line notation method that uses ASCII strings to represent
    molecules. More about SMILES can be found in the following link
    https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system

    **parameters:**
        obmol : openbabel molecule object

    **returns:**
        smi: SMILES notatation of the molecule.
    '''
    conv = ob.OBConversion()
    conv.SetOutFormat("smi")
    smi = conv.WriteString(obmol).rstrip()
    return smi


def ase_2_pybel(atoms):
    """
    As simple script to convert from ase atom object to pybel. There are
    many functionalities like atom typing, bond typing and autmatic addition of
    hydrogen that can be performed on a pybel molecule object, which can not be
    directly performed on an ase atom object.
    E.g\n
    add_hygrogen = pybel.addh()\n
    remove_hydrogen = pybel.removeh()\n
    https://openbabel.org/docs/dev/UseTheLibrary/Python_Pybel.html

    **parameters:**
        atoms : ASE atoms object

    **returns:**
        pybel: pybel molecular object.
    """
    a_str = ase_2_xyz(atoms)
    pybel_mol = pb.readstring("xyz", a_str)

    return pybel_mol


def max_index(lists):
    '''
    Extract index of list with the maximum element
    '''
    for i, elt in enumerate(lists):
        if elt == max(lists):
            return i


def compute_openbabel_cheminformatic(ase_atom):
    '''
    A function the returns all smiles, inchi and inchikeys
    from an ase_atom object. The function starts by wrapping
    the system, which is important for systems in which minimum
    image convertion has made atoms appear to be uncoordinated
    at random positions. The wrap functions coordinates all the
    atoms together.

    **parameters:**
        atoms : ASE atoms object

    **returns:**
        smi, inChi, inChiKey
    '''

    new_ase_atom = wrap_systems_in_unit_cell(ase_atom)

    new_ase_atom.set_pbc(False)
    pybel_mol = ase_2_pybel(new_ase_atom)

    obmol = pybel_mol.OBMol
    inChi,  inChiKey = compute_inchis(obmol)
    smi = compute_smi(obmol)
    return smi, inChi, inChiKey


def compute_cheminformatic_from_rdkit(ase_atom):
    '''
    A function that converts and ase atom into rkdit mol to extract some cheminformatic
    information.The script begins by taking an ase atom an creating an xyz string
    of the molecule, which is stored to memory. The xyz string is then read using
    the MolFromXYZBlock function in rdkit, which reads an xyz string of the molecule.
    However, this does not include any of the bonding information. To do this the following
    commands are parsed to the rdkit mol,bond_moll = Chem.Mol(rdmol)
    rdDetermineBonds.DetermineBonds(bond_moll,charge=0)
    https://github.com/rdkit/UGM_2022/blob/main/Notebooks/Landrum_WhatsNew.ipynb

    **parameters:**
        ASE atoms

    **returns:**
        smile string, inchi, inchikey
    '''
    if any(ase_atom.get_pbc()):
        ase_atom.set_pbc(False)

    xyz_string = ase_2_xyz(ase_atom)
    rdxyz = Chem.MolFromXYZBlock(xyz_string)
    rdmol = Chem.Mol(rdxyz)
    rdDetermineBonds.DetermineConnectivity(rdmol)
    rdDetermineBonds.DetermineBonds(rdmol, charge=0)
    smi = Chem.MolToSmiles(Chem.RemoveHs(rdmol))
    inChi = Chem.inchi.MolToInchi(rdmol)
    inChiKey = Chem.inchi.MolToInchiKey(rdmol)
    return smi, inChi, inChiKey


def compute_ase_neighbour(ase_atom):
    '''
    Create a connectivity graph using ASE neigbour list.

    **parameters:**
        ASE atoms

    **returns:**
        1.  atom_neighbors:
            A python dictionary, wherein each atom index
            is key and the value are the indices of it neigbours.
            e.g.
            atom_neighbors ={0:[1,2,3,4], 1:[3,4,5]...}

        2.  matrix:
            An adjacency matrix that wherein each row correspond to
            to an atom index and the colums correspond to the interaction
            between that atom to the other atoms. The entries in the
            matrix are 1 or 0. 1 implies bonded and 0 implies not bonded.

    '''
    atom_neighbors = {}
    cut_off = neighborlist.natural_cutoffs(ase_atom)

    neighbor_list = neighborlist.NeighborList(cut_off,
                                              self_interaction=False,
                                              bothways=True)
    neighbor_list.update(ase_atom)
    matrix = neighbor_list.get_connectivity_matrix(sparse=False)

    for atoms in ase_atom:
        connectivity, _ = neighbor_list.get_neighbors(atoms.index)
        atom_neighbors[atoms.index] = connectivity

    return atom_neighbors, matrix


def matrix2dict(bond_matrix):
    '''
    A simple procedure to convert an adjacency matrix into
    a python dictionary.

    **parameters:**
        bond matrix : adjacency matrix, type: nxn ndarray

    **returns:**
        graph: python dictionary
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

    **parameters:**
        graph: any python dictionary
        temp: a python list to hold nodes that have been visited
        node: a key in the python dictionary (graph), which is used as the starting or root node
        visited: python list containing nodes that have been traversed.

    **returns:**
        python dictionary
    '''
    visited[node] = True
    temp.append(node)
    for i in graph[node]:
        if visited[i] is False:
            temp = dfsutil_graph_method(graph, temp, i, visited)
    return temp


def longest_list(lst):
    '''
    return longest list in list of list
    '''
    return max(lst, key=len)


def remove_unbound_guest_and_return_unique(ase_atom):
    '''
    A simple script to remove guest from a metal organic framework.
    1)It begins by computing a connected graph component of all the fragments in the system
    using ASE neighbour list.
    2) Secondly it selects indicies of connected components which contain a metal
    3)if the there are two or more components, we create a pytmagen graph for each components and filter out all components that are not polymeric
    4) If there are two or more polymeric components, we check wether these systems there are identical or different and select only unique polymeric components

    **parameters:**
        ASE atoms

    **returns:**

        mof_indices : indinces of the guest free system. The guest free ase_atom object
        can be obtain as follows; E.g. guest_free_system = ase_atom[mof_indices]
    '''
    atom_neighbors, _ = compute_ase_neighbour(ase_atom)
    fragments = connected_components(atom_neighbors)
    if len(fragments) == 1:
        return [atom.index for atom in ase_atom]
    else:
        polymeric_indices = []
        for i in range(len(fragments)):
            super_cell = ase_atom[fragments[i]] * (2, 1, 1)
            coordination_graph, _ = compute_ase_neighbour(super_cell)
            pymat_graph = connected_components(coordination_graph)
            if len(pymat_graph) == 1:
                polymeric_indices.append(i)
        if len(polymeric_indices) > 0:
            Graphs = [StructureGraph.from_local_env_strategy(AseAtomsAdaptor.get_structure(
                ase_atom[fragments[i]]), JmolNN()) for i in polymeric_indices]
            temp_indices = [polymeric_indices[0]]
            unique = [Graphs[0]]
            for k in range(len(polymeric_indices[1:])):
                graph = Graphs[k]
                if True in [graph.diff(j) for j in unique]:
                    unique.append(Graphs[k])
                    temp_indices.append(polymeric_indices[i])
            mof_indices = []
            for frag_indices in temp_indices:
                mof_indices.extend(fragments[frag_indices])
            if len(mof_indices) == 0:
                return longest_list(fragments)
            else:
                return mof_indices
        else:
            return sum(fragments, [])


def remove_unbound_guest(ase_atom):
    '''
    A simple script to remove guest from a metal organic framework.
    1) It begins by computing a connected graph component of all the fragments in the system using ASE neighbour list.
    2) Secondly it selects indices of connected components which contain a metal
    3) if the there are two or more components, we create a pytmagen graph for each components and filter out all components that are not polymeric
    4) If there are two or more polymeric components, we check whether these systems there are identical or different and select all polymeric components

    **parameters:**
        ASE atoms

    **returns:**
        mof_indices: indices of the guest-free system. The guest-free ase_atom object
        can be obtained as follows; E.g.
        guest_free_system = ase_atom[mof_indices]
    '''
    atom_neighbors, _ = compute_ase_neighbour(ase_atom)
    fragments = connected_components(atom_neighbors)
    if len(fragments) == 1:
        return [atom.index for atom in ase_atom]
    else:
        polymeric_indices = []
        for i in range(len(fragments)):
            super_cell = ase_atom[fragments[i]] * (2, 1, 1)
            coordination_graph, _ = compute_ase_neighbour(super_cell)
            pymat_graph = connected_components(coordination_graph)
            if len(pymat_graph) == 1:
                polymeric_indices.append(i)
        if len(polymeric_indices) > 0:
            mof_indices = []
            for poly_index in polymeric_indices:
                mof_indices.extend(fragments[poly_index])
            if len(mof_indices) == 0:
                return longest_list(fragments)
            else:
                return mof_indices
        else:
            return sum(fragments, [])


def connected_components_recursive(graph):
    '''
    Find the connected fragments in a graph. Should work for any graph defined as a dictionary

    **parameters:**
        A graph in the form of dictionary e.g.::
        graph = {1:[0,1,3], 2:[2,4,5]}

    **returns:**
        returns a python list of list of connected components
        These correspond to individual molecular fragments.
        list_of_connected_components = [[1,2],[1,3,4]]
    '''
    visited = []
    list_of_connected_components = []
    for _ in list(graph.keys()):
        visited.append(False)
    for v in list(graph.keys()):
        if visited[v] is False:
            temp = []
            list_of_connected_components.append(
                dfsutil_graph_method(graph, temp, v, visited))
    return list_of_connected_components


def dfs_iterative(graph, start, visited):
    '''
    Depth-first search (DFS) iterative graph algorithm for traversing graph data structures.
    It starts at the root node 'start' and explores as far as possible along each branch before backtracking.
    This function is used as a utility for identifying connected components in the MOF graph.

    **parameters:**
        graph: a python dictionary where keys represent nodes and values are lists of neighboring nodes.
        start: a key in the python dictionary (graph), which is used as the starting or root node for the DFS.
        visited: a dictionary containing nodes that have been visited (True if visited, False otherwise).

    **returns:**
        A python list containing the nodes in the connected component starting from 'start'.
    '''
    component = []
    stack = [start]

    while stack:
        node = stack.pop()
        if not visited[node]:
            visited[node] = True
            component.append(node)
            # Add all unvisited neighbors to the stack.
            for neighbor in graph.get(node, []):
                if not visited[neighbor]:
                    stack.append(neighbor)

    return component


def connected_components_iterative(graph):
    '''
    Find the connected fragments in a graph. Should work for any graph defined as a dictionary.

    **parameters:**
        graph: a python dictionary representing the graph where keys are nodes and values are lists of adjacent nodes.

    **returns:**
        A python list of lists containing the connected components.
        Each sublist corresponds to an individual connected component (i.e., a set of nodes that are all connected).
    '''
    visited = {node: False for node in graph}
    components = []

    for node in graph:
        if not visited[node]:
            comp = dfs_iterative(graph, node, visited)
            components.append(comp)

    return components

def connected_components(graph):
    '''
    Selects the appropriate connected components algorithm based on the size of the system.
    Generally, the recursive method is faster but may exceed the recursion depth limit when the system is too large,
    while the iterative approach is slower but prevents the code from crashing in such cases.
    This function chooses the connected components method based on the size of the graph.

    **parameters:**
        graph: A dictionary representing a graph, where keys are nodes and values are lists of adjacent nodes.

    **returns:**
        A list of lists, where each sublist represents a connected component in the graph.
    '''
    if len(graph) > 1000:
        return connected_components_iterative(graph)
    else:
        return connected_components_recursive(graph)

def check_planarity(p1, p2, p3, p4):
    '''
    A simple procedure to check whether a point is planar to three other points.
    Important to distinguish porphyrin type metals

    **parameters:**
        p1, p2, p3, p4 : ndarray containing x,y,z values

    **returns:**
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
    # \------------------------------------------------
    a = b1 * c2 - b2 * c1
    b = a2 * c1 - a1 * c2
    c = a1 * b2 - b1 * a2
    d = round((-a * p1[0] - b * p1[1] - c * p1[2]), 0)
    # \------------------------------------------------
    factor = round((a * p4[0] + b * p4[1] + c * p4[2]), 0)
    verify = factor + d
    # \------------------------------------------------
    if verify == 0:
        planar = True
    return planar


def metal_in_porphyrin(ase_atom, graph):
    '''
    A funtion to check whether a metal is found at the centre of a phorphirin.
    The function specifically identifies metal atoms that
    are coordinated with four nitrogen atoms, as typically seen in porphyrin complexes.
    These atoms must also satisfy the planarity condition of the porphyrin structure.

    https://en.wikipedia.org/wiki/Transition_metal_porphyrin_complexes

    **parameters:**
        ase_atom: ASE atom
        graph: python dictionary containing neigbours

    **returns:**
        list of indices consisting of index of metal atoms found in the ASE atom
    '''
    all_porphyrin = []
    all_metal_symbols = [
        atom.index for atom in ase_atom if atom.symbol in transition_metals()]
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

    **parameters:**
        ase_atom: ASE atom
        graph: python dictionary containing neigbours

    **returns:**
        list of indices consisting of index of metal atoms found in the ASE atom
    '''
    old_list_of_connected_components = connected_components(graph)
    all_porphyrin = []
    all_metal_symbols = [
        atom.index for atom in ase_atom if atom.symbol in transition_metals()]
    metal_tmp = []
    N_tmp = []
    for idx in all_metal_symbols:
        connected = graph[idx]
        all_nitrogens = [i for i in connected if ase_atom[i].symbol == 'N']
        if len(all_nitrogens) == 4:
            metal_tmp.append(idx)
            N_tmp.extend(all_nitrogens)
    new_atom_indices = [i.index for i in ase_atom if i.index not in metal_tmp]
    tmp_atom = ase_atom[new_atom_indices]
    atom_neighbors, _ = compute_ase_neighbour(tmp_atom)
    list_of_connected_components = connected_components(atom_neighbors)
    if len(list_of_connected_components) == len(old_list_of_connected_components):
        all_porphyrin = metal_tmp + N_tmp
    return all_porphyrin


def move2front(index_value, coords):
    '''
    Move an index from any position in the list to the front
    The function is important to set the cell of a rodmof to point in thea-axis. Such that the system can be grow along this axis

    **parameters:**
        - index_value: index of item to move to the from
        - coords: list of coords

    **returns:**
        list of coords
    '''
    if any(isinstance(el, list) for el in coords):
        for data in coords:
            data.insert(0, data.pop(index_value))
    else:
        coords.insert(0, coords.pop(index_value))
    return coords


def find_carboxylates(ase_atom, graph):
    '''
    A simple aglorimth to search for carboxylates found in the system.

    **parameters:**
        - ase_atom: ASE atom

    **returns:**
        - dictionary of key = carbon index and values = oxygen index
    '''
    carboxyl = {}
    for atoms in ase_atom:
        if atoms.symbol == 'C':
            index = atoms.index
            oxygen = [i for i in graph[index] if ase_atom[i].symbol == 'O']
            if len(oxygen) == 2:
                oxy_metal = sum(
                    [[j for j in graph[i] if ase_atom[j].symbol in transition_metals()] for i in oxygen], [])
                if len(oxy_metal) > 0:
                    carboxyl[index] = oxygen
    return carboxyl


def find_carbonyl_sulphate(ase_atom, graph):
    """
    A simple aglorimth to search for Carbonyl sulphate  found in the system.
    ::

        O
        |
     -C-S
        |
        O

    **parameters:**
        ase_atom: ASE atom

    **returns:**
        dictionary of key = carbon index and values = oxygen index
    """
    sulphate = {}
    for atoms in ase_atom:
        if atoms.symbol == 'S':
            index = atoms.index
            oxygen = [i for i in graph[index] if ase_atom[i].symbol == 'O']
            carbon = [i for i in graph[index] if ase_atom[i].symbol == 'C']
            if len(oxygen) >= 1:
                oxy_metal = sum(
                    [[j for j in graph[i] if ase_atom[j].symbol in transition_metals()] for i in oxygen], [])
                if len(oxy_metal) > 0:
                    if len(carbon) > 0:
                        sulphate[index] = oxygen
    return sulphate


def find_sulfides(ase_atom, graph):
    """
    A simple aglorimth to search for sulfides.
    ::

        S
        |
       -C
        |
        S

    **parameters:**
        ase_atom: ASE atom

    **returns:**
        dictionary of key = carbon index and values = sulphur index
    """
    sulfides = {}
    metals = {}
    for atoms in ase_atom:
        if atoms.symbol == 'C':
            index = atoms.index
            sulphure_atoms = [i for i in graph[index]
                              if ase_atom[i].symbol == 'S']
            if len(sulphure_atoms) == 2:
                sulphure_to_metal = sum(
                    [[j for j in graph[i] if ase_atom[j].symbol in transition_metals()] for i in sulphure_atoms], [])
                if len(sulphure_to_metal) > 0:
                    sulfides[index] = sulphure_atoms
                    metals[index] = sulphure_to_metal
    return sulfides


def find_phosphite(ase_atom, graph):
    '''
    A simple aglorimth to search for sulfides.
    ::

        P
        |
       -C
        |
        P

    **parameters:**
        ase_atom: ASE atom

    **returns:**
        dictionary of key = carbon index and values = phosphorous index
    '''
    phosphorous = {}
    for atoms in ase_atom:
        if atoms.symbol == 'C':
            index = atoms.index
            phosphorous_atoms = [i for i in graph[index]
                                 if ase_atom[i].symbol == 'P']
            if len(phosphorous_atoms) == 2:
                phosphorous_to_metal = sum(
                    [[j for j in graph[i] if ase_atom[j].symbol in transition_metals()] for i in phosphorous_atoms], [])
                if len(phosphorous_to_metal) > 0:
                    phosphorous[index] = phosphorous_atoms
    return phosphorous


def find_COS(ase_atom, graph):
    '''
    A simple aglorimth to search for COS.
    ::

        O
        |
        -C
        |
        S

    **parameters:**
        ase_atom: ASE atom

    **returns:**
        dictionary of key = carbon index and values = sulphur index
    '''
    sulfides = {}
    for atoms in ase_atom:
        if atoms.symbol == 'C':
            index = atoms.index
            if len(graph[index]) == 3:
                sulphure_atoms = [i for i in graph[index]
                                  if ase_atom[i].symbol == 'S']
                oxygen_atoms = [i for i in graph[index]
                                if ase_atom[i].symbol == 'O']
                if len(sulphure_atoms) + len(oxygen_atoms) == 2:
                    sulphure_to_metal = sum(
                        [[j for j in graph[i] if ase_atom[j].symbol in transition_metals()] for i in sulphure_atoms], [])
                    oxygen_to_metal = sum(
                        [[j for j in graph[i] if ase_atom[j].symbol in transition_metals()] for i in oxygen_atoms], [])
                    if len(sulphure_to_metal) > 0 and len(oxygen_to_metal) > 0:
                        sulfides[index] = sulphure_atoms + oxygen_atoms
    return sulfides


def find_phosphate(ase_atom, graph):
    '''
    A simple algorithm to search for Carbonyl sulphate found in the system.
    ::

       O
       |
      -P-o
       |
       O

    **parameters:**
        ase_atom: ASE atom

    **returns**
        dictionary of key = carbon index and values = oxygen index
    '''
    phosphate = {}
    for atoms in ase_atom:
        if atoms.symbol == 'P':
            index = atoms.index
            oxygen = [i for i in graph[index] if ase_atom[i].symbol == 'O']
            if len(oxygen) >= 1:
                oxy_metal = sum(
                    [[j for j in graph[i] if ase_atom[j].symbol in transition_metals()] for i in oxygen], [])
                if len(oxy_metal) > 0:
                    phosphate[index] = oxygen
    return phosphate


def secondary_building_units(ase_atom):
    """
    1) Search for all carboxylate that are connected to a metal.
       Cut at the position between the carboxylate carbon and the adjecent carbon.
    2) Find all Nitrogen connected to metal. Check whether the nitrogen is in the
       centre of a porphirin ring. If no, cut at nitrogen metal bond.
    3) Look for oxygen that is connected to metal and two carbon. cut at metal oxygen bond

    **parameters:**
        ase_atom: ASE atom

    **returns:**
        list_of_connected_components: list of connected components,in which each list contains atom indices
        atom_pairs_at_breaking_point : Dictionary containing point of disconnection Porphyrin_checker: Boolean showing whether the metal is in the centre of a porpherin
        Regions: Dictionary of regions.
    """
    graph, bond_matrix = compute_ase_neighbour(ase_atom)
    porphyrin_checker = metal_in_porphyrin2(ase_atom, graph)
    atom_pairs_at_breaking_point = {}
    all_regions = {}
    bonds_to_break = []
    carboxylates = find_carboxylates(ase_atom, graph)
    all_sulphates = find_carbonyl_sulphate(ase_atom, graph)
    all_phosphites = find_phosphite(ase_atom, graph)
    all_sulphites = find_sulfides(ase_atom, graph)
    cos_group = find_COS(ase_atom, graph)
    ferocene_metal = all_ferrocene_metals(ase_atom, graph)
    all_metals = [i.index for i in ase_atom if i.symbol in transition_metals()]
    all_metals = [
        i for i in all_metals if i not in ferocene_metal+porphyrin_checker]
    for atoms in graph:
        if atoms in list(carboxylates.keys()):
            connected = graph[atoms]
            all_carbon_indices = [
                i for i in connected if ase_atom[i].symbol == 'C']
            all_nitrogens = [i for i in connected if ase_atom[i].symbol == 'N']
            S_indx = [i for i in connected if ase_atom[i].symbol == 'S']
            if len(all_carbon_indices) == 1:
                bonds_to_break.append([atoms] + all_carbon_indices)
                atom_pairs_at_breaking_point[atoms] = all_carbon_indices[0]
            if len(all_nitrogens) == 1:
                bonds_to_break.append([atoms] + all_nitrogens)
                atom_pairs_at_breaking_point[atoms] = all_nitrogens[0]
            if len(S_indx) == 1:
                bonds_to_break.append([atoms] + S_indx)
                atom_pairs_at_breaking_point[atoms] = S_indx[0]

        if ase_atom[atoms].symbol == 'C':
            connected = graph[atoms]
            oxygens = [i for i in connected if ase_atom[i].symbol == 'O']
            carbon_metal = [
                i for i in connected if ase_atom[i].symbol in transition_metals()]
            carbon_metal = [i for i in carbon_metal if i not in ferocene_metal]
            if len(oxygens) == 1:
                oxy_metal = [
                    i for i in oxygens if ase_atom[i].symbol in transition_metals()]
                oxy_metal = [i for i in oxy_metal if i not in ferocene_metal]
                if len(oxy_metal) == 1:
                    atom_pairs_at_breaking_point[oxygens[0]] = oxy_metal[0]
                    bonds_to_break.append([oxygens] + oxy_metal)
            if len(carbon_metal) > 0:
                for met in carbon_metal:
                    atom_pairs_at_breaking_point[atoms] = met
                    bonds_to_break.append([atoms] + [met])

        if atoms in list(all_sulphates.keys()):
            connected = graph[atoms]
            all_carbon_indices = [
                i for i in connected if ase_atom[i].symbol == 'C']
            oxygen = all_sulphates[atoms]
            if len(all_carbon_indices) == 1:
                atom_pairs_at_breaking_point[atoms] = all_carbon_indices[0]
                bonds_to_break.append([atoms] + all_carbon_indices)
            if len(all_carbon_indices) > 1:
                for oxy in oxygen:
                    metal = [i for i in graph[oxy]
                             if ase_atom[i].symbol in transition_metals()]
                    metal = [i for i in metal if i not in ferocene_metal]
                    if len(metal) > 0:
                        for met in metal:
                            atom_pairs_at_breaking_point[oxy] = met
                            bonds_to_break.append([oxy] + [met])

        if atoms in list(all_phosphites.keys()):
            all__n_indices = all_phosphites[atoms]
            connected = [i for i in graph[atoms] if i not in all__n_indices]
            for neigbour in connected:
                atom_pairs_at_breaking_point[atoms] = neigbour
                bonds_to_break.append([atoms] + [neigbour])

        if atoms in list(all_sulphites.keys()):
            all__n_indices = all_sulphites[atoms]
            connected = [i for i in graph[atoms] if i not in all__n_indices]
            for neigbour in connected:
                atom_pairs_at_breaking_point[atoms] = neigbour
                bonds_to_break.append([atoms] + [neigbour])

        # if atoms in list(all_sulphites.keys()):
        #     all__n_indices = all_sulphites[atoms]
        #     neigbours = sum([graph[i].tolist() for i in all__n_indices], [])
        #     metals = [
        #         i for i in neigbours if ase_atom[i].symbol in transition_metals()]
        #     metals = [i for i in metals if i not in ferocene_metal]
        #     if len(metals) > 0:
        #         for met in metals:
        #             atom_pairs_at_breaking_point[atoms] = met
        #             bonds_to_break.append([atoms] + [met])

        if atoms in list(cos_group.keys()):
            connected = graph[atoms]
            all_carbon_indices = [
                i for i in connected if ase_atom[i].symbol == 'C']
            non_metals = cos_group[atoms]
            if len(all_carbon_indices) == 1:
                atom_pairs_at_breaking_point[atoms] = all_carbon_indices[0]
                bonds_to_break.append([atoms] + all_carbon_indices)
            if len(all_carbon_indices) > 1:
                for atom_idx in non_metals:
                    metal = [i for i in graph[atom_idx]
                             if ase_atom[i].symbol in transition_metals()]
                    metal = [i for i in metal if i not in ferocene_metal]
                    if len(metal) > 0:
                        for met in metal:
                            atom_pairs_at_breaking_point[atom_idx] = met
                            bonds_to_break.append([atom_idx] + [met])

        if ase_atom[atoms].symbol == 'O':
            seen = sum(list(carboxylates.values())
                       + list(all_sulphates.values())
                       + list(cos_group.values()), [])
            if atoms not in seen:
                connected = graph[atoms]
                metal = [
                    i for i in connected if ase_atom[i].symbol in transition_metals()]
                metal = [i for i in metal if i not in ferocene_metal]
                Nitrogen = [i for i in connected if ase_atom[i].symbol == 'N']
                carbon = [i for i in connected if ase_atom[i].symbol
                          == 'C' and i not in list(carboxylates.keys())]
                if len(metal) >= 1 and len(carbon) == 1:
                    atom_pairs_at_breaking_point[atoms] = carbon[0]
                    bonds_to_break.append([atoms] + carbon)
                if len(metal) == 1 and len(carbon) == 2:
                    atom_pairs_at_breaking_point[atoms] = metal[0]
                    bonds_to_break.append([atoms] + metal)
                if len(metal) == 1 and len(Nitrogen) == 1:
                    n_carbon = [i for i in graph[Nitrogen[0]]
                                if ase_atom[i].symbol == 'C' and i not in list(carboxylates.keys())]
                    n_nitrogen = [i for i in graph[Nitrogen[0]]
                                  if ase_atom[i].symbol == 'N']
                    n_sulphur = [i for i in graph[Nitrogen[0]]
                                 if ase_atom[i].symbol == 'S' and i not in list(all_sulphates.keys())]
                    if len(n_carbon) > 1:
                        atom_pairs_at_breaking_point[atoms] = Nitrogen[0]
                        bonds_to_break.append([atoms] + Nitrogen)

                    elif len(n_nitrogen) > 1:
                        atom_pairs_at_breaking_point[atoms] = Nitrogen[0]
                        bonds_to_break.append([atoms] + Nitrogen)
                    elif len(n_sulphur) > 1:
                        atom_pairs_at_breaking_point[atoms] = Nitrogen[0]
                        bonds_to_break.append([atoms] + Nitrogen)

        if ase_atom[atoms].symbol == 'N':
            connected = graph[atoms]
            metal = [
                i for i in connected if ase_atom[i].symbol in transition_metals()]
            metal = [i for i in metal if i not in ferocene_metal]
            if atoms not in porphyrin_checker:
                if len(metal) == 1:
                    atom_pairs_at_breaking_point[atoms] = metal[0]
                    bonds_to_break.append([atoms] + metal)

        if ase_atom[atoms].symbol == 'S':
            seen = sum(list(carboxylates.values())
                       + list(all_sulphates.values())
                       + list(all_sulphites.values())
                       + list(cos_group.values()), []
                       )
            if atoms not in seen:
                connected = graph[atoms]
                metal = [
                    i for i in connected if ase_atom[i].symbol in transition_metals()]
                metal = [i for i in metal if i not in ferocene_metal]
                if len(metal) > 0:
                    for met in metal:
                        atom_pairs_at_breaking_point[atoms] = met
                        bonds_to_break.append([atoms, met])

        if ase_atom[atoms].symbol == 'P':
            # Find the carbon closest to P, which is not bonded to a metal and cut
            # 1) Look for phosphorous atoms
            # 2) Look for all it's neigbours
            # 3) Look for neigbours that are not connected to metal or hydrogen.
            seen = sum(list(all_phosphites.values()), [])
            if atoms not in seen:
                connected = graph[atoms]
                # not_connected_to_metal_or_hygrogen = [[i for i in graph[j] if ase_atom[i].symbol not in transition_metals() or ase_atom[i].symbol != 'H'] for j in connected]

                metal_oxy = [[i for i in graph[j] if ase_atom[i].symbol in transition_metals()]
                             for j in connected]

                metal = sum(metal_oxy, [])
                metal = [i for i in metal if i not in ferocene_metal]
                closest_atoms = sum(
                    [[i for i in graph[j] if i != atoms and not ase_atom[i].symbol in transition_metals()] for j in connected], [])

                if len(metal) > 0:
                    all_carbon_indices = sum([[i for i in graph[j] if i in connected]
                                              for j in closest_atoms], [])
                    for frag in all_carbon_indices:
                        atom_pairs_at_breaking_point[atoms] = frag
                        atom_pairs_at_breaking_point[frag] = atoms
                    bonds_to_break.append([atoms, frag])

        if ase_atom[atoms].symbol == 'B':
            # Find the carbon closest to P, which is  bonded to a metal and cut
            connected = [i for i in graph[atoms]
                         if ase_atom[i].symbol not in transition_metals()]
            metal_oxy = [[[i, j] for i in graph[j] if ase_atom[i].symbol in transition_metals()]
                         for j in connected]

            metal_connect = sum(metal_oxy, [])
            metal = [i for i in metal if i not in ferocene_metal]
            if len(metal_connect) > 0:
                for frag in metal_connect:
                    atom_pairs_at_breaking_point[frag[0]] = frag[1]
                    bonds_to_break.append(frag)
    # In special cases some carbon and hydrogen get very closed to metals
    # So in such case it is important to clear this
    for metal in all_metals:
        connected = graph[metal]
        c_H = [i for i in connected if ase_atom[i].symbol in ['H', 'C']]
        if len(c_H) > 0:
            for c_h in c_H:
                bonds_to_break.append([c_h, metal])
                atom_pairs_at_breaking_point[metal] = c_h

    for bonds in bonds_to_break:
        bond_matrix[bonds[0], bonds[1]] = 0
        bond_matrix[bonds[1], bonds[0]] = 0

    new_ase_graph = matrix2dict(bond_matrix)
    try:
        list_of_connected_components = connected_components(new_ase_graph)
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
        atom_pairs_at_breaking_point,
        porphyrin_checker, all_regions
    ]


def ligands_and_metal_clusters(ase_atom):
    '''
    Start by checking whether there are more than 2 layers
    if yes, select one
    Here we select the largest connected component

    **parameters:**
        ase_atom: ASE atom

    **returns:**
        list_of_connected_components: list of connected components, in which each list contains atom indices
        atom_pairs_at_breaking_point: Dictionary containing point of disconnection
        Porpyrin_checker: Boolean showing whether the metal is in the centre of a porpherin
        Regions: Dictionary of regions.
    '''
    graph, bond_matrix = compute_ase_neighbour(ase_atom)

    porphyrin_checker = metal_in_porphyrin2(ase_atom, graph)

    all_regions = {}
    atom_pairs_at_breaking_point = {}
    bonds_to_break = []
    seen_oxygen = []
    carboxylates = find_carboxylates(ase_atom, graph)
    all_sulphates = find_carbonyl_sulphate(ase_atom, graph)
    all_phosphites = find_phosphite(ase_atom, graph)
    all_sulphites = find_sulfides(ase_atom, graph)
    ferocene_metal = all_ferrocene_metals(ase_atom, graph)
    all_metals = [i.index for i in ase_atom if i.symbol in transition_metals()]
    all_metals = [
        i for i in all_metals if i not in ferocene_metal+porphyrin_checker]
    for atoms in graph:
        if atoms in list(carboxylates.keys()):
            oxygen = carboxylates[atoms]
            for oxy in oxygen:
                metal = [i for i in graph[oxy]
                         if ase_atom[i].symbol in transition_metals()]
                metal = [i for i in metal if i not in ferocene_metal]
                if len(metal) > 0:
                    seen_oxygen.append(oxy)
                    for met in metal:
                        bonds_to_break.append([atoms] + [met])
                        atom_pairs_at_breaking_point[atoms] = met

        if atoms in list(all_sulphates.keys()):
            oxygen = all_sulphates[atoms]
            connected = graph[atoms]
            # all_carbon_indices = [i for i in connected if ase_atom[i].symbol == 'C']
            for oxy in oxygen:
                metal = [i for i in graph[oxy]
                         if ase_atom[i].symbol in transition_metals()]
                metal = [i for i in metal if i not in ferocene_metal]
                if len(metal) > 0:
                    for met in metal:
                        bonds_to_break.append([oxy] + [met])
                        atom_pairs_at_breaking_point[oxy] = met

        if atoms in list(all_phosphites.keys()):
            all__n_indices = all_phosphites[atoms]
            for phos in all__n_indices:
                metals = [i for i in graph[phos]
                          if ase_atom[i].symbol in transition_metals()]
                for met in metals:
                    atom_pairs_at_breaking_point[phos] = met
                    bonds_to_break.append([phos, met])

        if atoms in list(all_sulphites.keys()):
            all__n_indices = all_sulphites[atoms]
            for sulphur in all__n_indices:
                metals = [i for i in graph[sulphur]
                          if ase_atom[i].symbol in transition_metals()]
                for met in metals:
                    atom_pairs_at_breaking_point[sulphur] = met
                    bonds_to_break.append([sulphur, met])

        if ase_atom[atoms].symbol == 'C':
            connected = graph[atoms]
            carbon_metal = [
                i for i in connected if ase_atom[i].symbol in transition_metals()]
            carbon_metal = [i for i in carbon_metal if i not in ferocene_metal]
            if len(carbon_metal) > 0:
                for met in carbon_metal:
                    atom_pairs_at_breaking_point[atoms] = met
                    bonds_to_break.append([atoms] + [met])

        if ase_atom[atoms].symbol == 'N':
            connected = graph[atoms]
            metal = [
                i for i in connected if ase_atom[i].symbol in transition_metals()]
            metal = [i for i in metal if i not in ferocene_metal]
            if len(metal) > 0 and atoms not in porphyrin_checker:
                for met in metal:
                    bonds_to_break.append([atoms, met])
                    atom_pairs_at_breaking_point[atoms] = met
                # atom_pairs_at_breaking_point [metal[0]] = atoms

        if ase_atom[atoms].symbol == 'S':
            seen = sum(list(all_sulphites.values()), [])
            if atoms not in seen:
                connected = graph[atoms]
                metal = [
                    i for i in connected if ase_atom[i].symbol in transition_metals()]
                metal = [i for i in metal if i not in ferocene_metal]
                if len(metal) > 0:
                    if len(metal) > 0:
                        for met in metal:
                            bonds_to_break.append([atoms] + [met])
                            atom_pairs_at_breaking_point[atoms] = met

        if ase_atom[atoms].symbol == 'O':
            # if not atoms in Seen_oxygen:
            connected = graph[atoms]
            metal = [
                i for i in connected if ase_atom[i].symbol in transition_metals()]
            metal = [i for i in metal if i not in ferocene_metal]
            Nitrogen = [i for i in connected if ase_atom[i].symbol == 'N']
            carbon = [i for i in connected if ase_atom[i].symbol == 'C']
            if len(metal) > 0 and len(carbon) == 1:
                for met in metal:
                    atom_pairs_at_breaking_point[atoms] = met
                    bonds_to_break.append([atoms, met])
            if len(metal) > 0 and len(Nitrogen) == 1:
                n_carbon = [i for i in graph[Nitrogen[0]]
                            if ase_atom[i].symbol in ['C', 'S', 'N']]
                if len(n_carbon) > 1:
                    for met in metal:
                        atom_pairs_at_breaking_point[atoms] = met
                        bonds_to_break.append([atoms, met])
            if len(carbon) > 1 and len(metal) > 0:
                for met in metal:
                    bonds_to_break.append([atoms] + [met])
                    atom_pairs_at_breaking_point[atoms] = met

        if ase_atom[atoms].symbol == 'P':
            # Find the carbon closest to P, which is not bonded to a metal and cut
            seen = sum(list(all_phosphites.values()), [])
            if atoms not in seen:
                connected = [i for i in graph[atoms]
                             if ase_atom[i].symbol not in transition_metals()]
                metal_oxy = [[i for i in graph[j] if ase_atom[i].symbol in transition_metals()]
                             for j in connected]
                metal = sum(metal_oxy, [])
                metal = [i for i in metal if i not in ferocene_metal]
                closest_atoms = sum([[[i, j] for i in graph[j] if i != atoms and ase_atom[i].symbol in
                                    transition_metals()]for j in connected], [])
                for frag in closest_atoms:
                    bonds_to_break.append(frag)
                    atom_pairs_at_breaking_point[frag[0]] = frag[1]
    # In special cases some carbon and hydrogen get very closed to metals
    # So in such case it is important to clear this
    for metal in all_metals:
        connected = graph[metal]
        c_H = [i for i in connected if ase_atom[i].symbol in ['H', 'C']]
        if len(c_H) > 0:
            for c_h in c_H:
                bonds_to_break.append([c_h, metal])
                atom_pairs_at_breaking_point[metal] = c_h

    for bonds in bonds_to_break:
        bond_matrix[bonds[0], bonds[1]] = 0
        bond_matrix[bonds[1], bonds[0]] = 0

    new_ase_graph = matrix2dict(bond_matrix)
    try:
        list_of_connected_components = connected_components(new_ase_graph)
    except Exception:
        import networkx as nx
        N_Graph = nx.from_dict_of_lists(new_ase_graph)
        list_of_connected_components = [
            list(i) for i in list(nx.connected_components(N_Graph))]

    all_pm_structures = [sorted(ase_atom[i].symbols)
                         for i in list_of_connected_components]
    for i in range(len(all_pm_structures)):
        temp = []
        for j in range(len(all_pm_structures)):
            if all_pm_structures[i] == all_pm_structures[j]:
                temp.append(j)
        if temp not in all_regions .values():
            all_regions[i] = temp

    return list_of_connected_components, atom_pairs_at_breaking_point, porphyrin_checker, all_regions


def is_rodlike(metal_sbu):
    '''
    Simple test to check whether a metal sbu is a rodlike MOF

    **parameter:**
        metal_sbu : ase_atom

    **returns:**
        bool : True if the metal sbu is a rodlike MOF, False otherwise
    '''
    rod_check = []
    cells = [(2, 1, 1), (1, 2, 1), (1, 1, 2)]
    for index, ijk in enumerate(cells):
        rod = metal_sbu * ijk
        graph, _ = compute_ase_neighbour(rod)
        list_of_connected_components = connected_components(graph)
        if len(list_of_connected_components) == 1:
            rod_check.append(index)
    return rod_check


def all_ferrocene_metals(ase_atom, graph):
    '''
    A function to find metals corresponding to ferrocene.
    These metals should not be considered during mof-constructions

    **parameters:**
        ase_atom: ASE atom
        graph : dictionary containing neigbour lists

    **returns:**
        list_of_metals: list of indices of ferrocene metals
    '''
    list_of_metals = []
    graph, _ = compute_ase_neighbour(ase_atom)
    for atom_index in graph:
        if ase_atom[atom_index].symbol in transition_metals():
            connectivity = graph[atom_index]
            if len(connectivity) >= 10:
                number_of_carbons = 0
                for neigbour in connectivity:
                    if ase_atom[neigbour].symbol == 'C':
                        number_of_carbons += 1
                if number_of_carbons >= 10:
                    list_of_metals.append(atom_index)
    return list_of_metals


def is_ferrocene(metal_sbu, graph):
    '''
    A simple script to check whether a metal_sbu is ferrocene

    **parameter:**
        metal_sbu : ase_atom
        graph : dictionary containing neigbour lists

    **returns:**
        bool : True if the metal_sbu is ferrocene, False otherwise
    '''
    check = []
    verdict = False
    all_connectivity = list(graph.values())
    for connectivity in all_connectivity:
        if len(connectivity) >= 10:
            carbons = 0
            for bondedAtomIndex in connectivity:
                if metal_sbu[bondedAtomIndex].symbol == 'C':
                    carbons += 1
            if carbons >= 10:
                verdict = True
            check.append(verdict)
        correct = False
        if True in check:
            correct = True
    return correct


def is_paddlewheel(metal_sbu, graph):
    """
    returns True if the atom is part of a paddlewheel motif

    **parameter:**
        metal_sbu : ase_atom
        graph : dictionary containing neigbour lists

    **returns:**
        bool : True if the metal_sbu is part of a paddlewheel motif, False otherwise
    """
    check = []
    verdict = False
    all_connectivity = list(graph.values())
    for connectivity in all_connectivity:
        metalNeighbours = 0
        oxygenNeighbours = 0
        for bondedAtomIndex in connectivity:
            if metal_sbu[bondedAtomIndex].symbol in transition_metals():
                metalNeighbours += 1
            if metal_sbu[bondedAtomIndex].symbol == 'O':
                oxygenNeighbours += 1
        if metalNeighbours == 1 and oxygenNeighbours == 4 and len(connectivity) >=\
                5 and len(connectivity) <= 6:
            verdict = True
        check.append(verdict)
    correct = False
    if True in check:
        correct = True
    return correct


def is_paddlewheel_with_water(ase_atom, graph):
    """
    returns True if the atom is part of a paddle wheel with water motif

    **parameter:**
        ase_atom : ase_atom
        graph : dictionary containing neigbour lists

    **returns:**
        bool : True if the metal_sbu is part of a paddle wheel with water motif, False otherwise
    """
    check = []
    metal = []
    verdict = False
    for atoms in ase_atom:
        if atoms.symbol in transition_metals():
            index = atoms.index
            metal.append(index)
            connectivity = graph[index]
            if len([ase_atom[i].symbol for i in connectivity if
                    ase_atom[i].symbol == 'O']) == 5:
                verdict = True
                check.append(verdict)
    correct = False
    if True in check and len(metal) == 2:
        correct = True
    return correct


def is_uio66(ase_atom, graph):
    """
    returns True if the atom is part of a UIO66 motif

    **parameter:**
        ase_atom : ase_atom
        graph : dictionary containing neigbour lists

    **returns:**
        bool : True if the metal_sbu is part of a UIO66 motif, False otherwise
    """
    check = []
    verdict = False
    all_connectivity = list(graph.values())
    for connectivity in all_connectivity:
        metal_neighbours = 0
        oxygen_neighbours = 0
        for bonded_atom_index in connectivity:
            if ase_atom[bonded_atom_index].symbol in transition_metals():
                metal_neighbours += 1
            if ase_atom[bonded_atom_index].symbol == 'O':
                oxygen_neighbours += 1
        if metal_neighbours == 4 and (oxygen_neighbours == 6 or oxygen_neighbours == 8)\
                and (len(connectivity) == 10 or len(connectivity) == 12):
            verdict = True
        check.append(verdict)
        correct = False
        if True in check:
            correct = True
    return correct


def is_irmof(ase_atom, graph):
    """
    returns True if the atom is part of a IRMOF motif

    **parameter:**
        ase_atom : ase_atom
        graph : dictionary containing neigbour lists

    **returns:**
        bool : True if the metal_sbu is part of a IRMOF motif, False otherwise
    """
    check = []
    verdict = False
    for atoms in ase_atom:
        if atoms.symbol == 'O':
            index = atoms.index
            connectivity = graph[index]
            if len(connectivity) == 4 and len([ase_atom[i].symbol for i in connectivity if
                                               ase_atom[i].symbol in transition_metals()]) == 4:

                verdict = True
                check.append(verdict)
    correct = False
    if True in check:
        correct = True
    return correct


def is_mof32(ase_atom, graph):
    """
    returns True if the atom is part of a MOF32 motif

    **parameter:**
        ase_atom : ase_atom
        graph : dictionary containing neigbour lists

    **returns:**
        bool : True if the metal_sbu is part of a MOF32 motif, False otherwise
    """
    check = []
    metal = []
    verdict = False
    for atoms in ase_atom:
        if atoms.symbol in transition_metals():
            index = atoms.index
            metal.append(index)
            connectivity = graph[index]
            if len([ase_atom[i].symbol for i in connectivity if ase_atom[i].symbol == 'O']) == 8:
                verdict = True
                check.append(verdict)
    correct = False
    if True in check and len(metal) == 1:
        correct = True
    return correct


def rod_manipulation(ase_atom, checker):
    '''
    Script to adjust Rodlike sbus.
    1) Its collects the axis responsible for expanding the rod
    2) It shifts all coordinates to the axis
    3) It rotates the rod to lie in the directions of expansion,

    **parameter:**
        ase_atom : ASE atoms object
        checker : list of indices of atoms responsible for rod expansion

    **returns:**
        ase_atom : ASE atoms object with adjusted coordinates
    '''

    cell = ase_atom.get_cell().tolist()
    value_index = checker[0]
    cell_value = cell[value_index]
    new_cell = move2front(checker, cell_value)
    coords = ase_atom.positions - cell_value
    new_position = move2front(checker, coords.tolist())
    ase_atom.positions = new_position
    return ase_atom, new_cell


def find_unique_building_units(list_of_connected_components, atom_pairs_at_breaking_point,
                               ase_atom, porphyrin_checker, all_regions, wrap_system=True, cheminfo=False, add_dummy=False):
    '''
    A function that identifies and processes unique building units within a MOF.
    This function deconstructs a MOF into its constituent building unit and identifies
    unique building units. It can also
    1. wrap the system into the unit cell,
    2. add dummy atoms to neutralize the building units,
    3. compute cheminformatics data using Open Babel.

    **parameters:**
        list_of_connected_components : list of connected components
        atom_pairs_at_breaking_point : dictionary of atom pairs at breaking point
        ase_atom : ASE atoms object
        porphyrin_checker : list of indices of porphyrins
        all_regions : dictionary of regions
        wrap_system : boolean, default True
        cheminfo : boolean, default False
        add_dummy: boolean, default False
        add_dummy keyword enables the addition of dumnmy atoms which can then be replaced
        by hydrogen to neutralize the building blocks.

    **returns:**
        mof_metal: list of metal indices
        mof_linker: list of linker indices
        concentration: dictionary of concentration of building units
    '''
    mof_metal = []
    mof_linker = []
    # concentration = {}
    building_unit_regions = {}
    for idx, key in enumerate(all_regions.keys()):
        frag = list(all_regions .keys())[idx]
        components = list_of_connected_components[all_regions[frag][0]]
        all_breaking_point = list(atom_pairs_at_breaking_point.keys(
        )) + list(atom_pairs_at_breaking_point.values())
        point_of_extension = [i for i in all_breaking_point if i in components]
        mapped_indices = dict(
            [(i, j) for i, j in zip(components, range(len(components)))])
        molecule_to_write = ase_atom[components]
        molecule_to_write.info['point_of_extension'] = [
            mapped_indices[i] for i in point_of_extension]

        molecule_for_prop = molecule_to_write
        if wrap_system and add_dummy:
            dummy_idx = [find_key_or_value(
                i, atom_pairs_at_breaking_point) for i in point_of_extension]

            if dummy_idx:
                dummy_mol = ase_atom[dummy_idx]
                len_mole_without_dummy = len(molecule_to_write)

                molecule_to_write += dummy_mol
                molecule_to_write = wrap_systems_in_unit_cell(
                    molecule_to_write)

                for i in range(len_mole_without_dummy, len(molecule_to_write)):
                    molecule_to_write[i].symbol = 'X'

                molecule_for_prop = molecule_to_write[:len_mole_without_dummy]

        if wrap_system is True and add_dummy is False:
            molecule_to_write = wrap_systems_in_unit_cell(molecule_to_write)

        if cheminfo:
            smi, chem_inchi, chem_inchiKey = compute_openbabel_cheminformatic(
                molecule_for_prop)
            molecule_to_write.info['smi'] = smi
            molecule_to_write.info['inchi'] = str(chem_inchi)
            molecule_to_write.info['inchikey'] = str(chem_inchiKey)

            molecule_to_write.info['inchikey'] = str(chem_inchiKey)
        atom_indices_mapping = [
            list_of_connected_components[i] for i in all_regions[key]]
        molecule_to_write.info['atom_indices_mapping'] = atom_indices_mapping

        # if add_dummy:
        #     dummy_idx = [find_key_or_value(i, atom_pairs_at_breaking_point) for i in  point_of_extension]
        #     dummy_mol = ase_atom[dummy_idx]
        #     # if wrap_system:
        #     #     dummy_mol = wrap_systems_in_unit_cell(dummy_mol)
        #     for atom in dummy_mol:
        #         atom.symbol = 'X'
        #     molecule_to_write += dummy_mol

        metal = [i.index for i in molecule_for_prop if i.symbol in transition_metals(
        ) and i.index not in porphyrin_checker]
        non_ferocene_metal = []
        if len(metal) > 0:
            graph_sbu, _ = compute_ase_neighbour(molecule_for_prop)
            # Check whether the metal sbu is a rod mof. If is it is rod mof,
            # we rotate and aligne the sbu such that the axis of rotation will be the a-axis.
            if len(is_rodlike(molecule_for_prop)) == 1:
                molecule_to_write.info['sbu_type'] = 'rodlike'
            elif is_paddlewheel(molecule_for_prop, graph_sbu):
                molecule_to_write.info['sbu_type'] = 'paddlewheel'
            elif is_paddlewheel_with_water(molecule_for_prop, graph_sbu):
                molecule_to_write.info['sbu_type'] = 'paddlewheel_with_water'
            elif is_uio66(molecule_for_prop, graph_sbu):
                molecule_to_write.info['sbu_type'] = 'UIO66_sbu'
            elif is_mof32(molecule_for_prop, graph_sbu):
                molecule_to_write.info['sbu_type'] = 'MOF32_sbu'
            elif is_irmof(molecule_for_prop, graph_sbu):
                molecule_to_write.info['sbu_type'] = 'IRMOF_sbu'
            elif is_ferrocene(molecule_for_prop, graph_sbu):
                molecule_to_write.info['sbu_type'] = 'ferrocenelike'
                ferocene_metal = all_ferrocene_metals(
                    molecule_for_prop, graph_sbu)
                non_ferocene_metal = [
                    i for i in metal if not i in ferocene_metal]
                if len(non_ferocene_metal) == 0:
                    mof_linker.append(molecule_to_write)
                    continue
            else:
                molecule_to_write.info['sbu_type'] = 'still checking!'
            mof_metal.append(molecule_to_write)

        else:
            mof_linker.append(molecule_to_write)
        building_unit_regions[idx] = molecule_to_write.info['atom_indices_mapping']
    return mof_metal, mof_linker, building_unit_regions


def metal_coordination_number(ase_atom):
    '''
    Extract coordination number of central metal

    **paramters:**
        ase_atom : ASE atoms object

    **returns:**
        metal_elt : list of metal elements
    '''
    metal_coordination = {}
    graph, _ = compute_ase_neighbour(ase_atom)
    porph_indices = metal_in_porphyrin2(ase_atom, graph)
    metal_indices = [
        i.index for i in ase_atom if i.symbol in transition_metals()]
    metal_elt = []
    for i in metal_indices:
        if i not in porph_indices and ase_atom[i].symbol not in metal_elt:
            metal_elt.append(ase_atom[i].symbol)
            metal_coordination[ase_atom[i].symbol] = len(graph[i])
    return metal_elt, metal_coordination


def metal_coordination_enviroment(ase_atom):
    '''
    Find the enviroment of a metal

    **paramters:**
        ase_atom : ASE atoms object

    **returns:**
        metal_elt : list of metal elements
        metal_enviroment : dict where key is metal element and value is list of neighboring elements
    '''
    metal_enviroment = {}
    seen = []
    graph, _ = compute_ase_neighbour(ase_atom)
    porph_indices = metal_in_porphyrin2(ase_atom, graph)
    metal_indices = [
        i.index for i in ase_atom if i.symbol in transition_metals()]
    for i in metal_indices:
        if i not in porph_indices and i not in seen:
            seen.append(ase_atom[i].symbol)
            metal_enviroment[ase_atom[i].symbol] = [
                ase_atom[j].symbol for j in graph[i]]
    return metal_enviroment


def mof_regions(ase_atom, list_of_connected_components, atom_pairs_at_breaking_point):
    '''
    A function to map all atom indices to exact position in which the find themselves in the MOF.
    This function is used to partition a MOF into regions that correspond to unique
    unique building units.

    **parameters:**
        ase_atom: ASE atom
        list_of_connected_components : list of list, wherein each list correspond to atom indices of a specific building unit
        atom_pairs_at_breaking_point: dictionary containing pairs of atoms from which the bonds were broken

    **returns:**
        Move an index from any position in the list to the front
        The function is important to set the cell of a rodmof to point in the
        a-axis. Such that the system can be grow along this axis
    '''
    all_regions = {}
    xis_regions = {}
    all_pm_structures = [sorted(ase_atom[i].symbols)
                         for i in list_of_connected_components]
    for i in range(len(all_pm_structures)):
        temp = []
        for j in range(len(all_pm_structures)):
            if all_pm_structures[i] == all_pm_structures[j]:
                temp.append(j)
        if not temp in all_regions.values():
            all_regions[i] = temp

    for idx in range(len(all_regions .keys())):
        frag = list(all_regions .keys())[idx]
        components = list_of_connected_components[all_regions[frag][0]]
        dummy = [[atom_pairs_at_breaking_point[j] for j in list(atom_pairs_at_breaking_point.keys()) if j in comp]
                 for comp in components]
        xis_regions[idx] = dummy
    return all_regions, xis_regions


def wrap_systems_in_unit_cell(ase_atom, max_iter=30, skin=0.3):
    '''
    A simple aglorithm to reconnnect all atoms wrapped in a periodic
    boundary condition such that all atoms outside the box will appear reconnected.

    **parameters:**
        ase_atom : ASE atoms object
        max_iter : int, optional (default=30) Maximum number of iterations for reconnection
        skin : float, optional (default=0.3) Skin distance for bond reconnection

    **returns:**
        - ase_atom : ASE atoms object with reconnected atoms wrapped in a periodic boundary condition
    '''
    if not any(ase_atom.get_pbc()):
        return ase_atom
    else:
        new_position = geometry.wrap_positions(ase_atom.positions, ase_atom.cell, pbc=True, center=(
            0, 0, 0), pretty_translation=True, eps=1e-07)
        ase_atom.positions = new_position

        graph, bond_matrix = compute_ase_neighbour(ase_atom)

        for atom in graph:
            connected = graph[atom]
            for nl in connected:
                check = covalent_radius(
                    ase_atom[atom].symbol) + covalent_radius(ase_atom[nl].symbol) + skin
                bond = round(ase_atom.get_distance(atom, nl), 2)
                if bond > check:
                    bond_matrix[atom][nl] = 0

        new_ase_graph = matrix2dict(bond_matrix)
        list_of_connected_components = connected_components(new_ase_graph)
        number_of_iterations = 0
        while len(list_of_connected_components) != 1:
            all_len = [len(i) for i in list_of_connected_components]
            max_index = all_len.index(max(all_len))
            Root = list_of_connected_components[max_index]
            list_of_connected_components.pop(max_index)
            All_sum = sum(list_of_connected_components, [])
            for atom in Root:
                connected = graph[atom]
                for nl in All_sum:
                    if nl in connected:
                        v = ase_atom[nl].position - ase_atom[atom].position
                        mic_vector = geometry.find_mic(
                            v, ase_atom.get_cell(), pbc=True)
                        ase_atom[nl].position = mic_vector[0] + \
                            ase_atom[atom].position

            graph, bond_matrix = compute_ase_neighbour(ase_atom)
            for atom in graph:
                connected = graph[atom]
                for nl in connected:
                    check = covalent_radius(
                        ase_atom[atom].symbol) + covalent_radius(ase_atom[nl].symbol) + skin
                    bond = round(ase_atom.get_distance(atom, nl), 2)
                    if bond > check:
                        bond_matrix[atom][nl] = 0
            new_ase_graph = matrix2dict(bond_matrix)
            list_of_connected_components = connected_components(new_ase_graph)
            number_of_iterations += 1
            if number_of_iterations == max_iter:
                break
        return ase_atom


def angle_tolerance_to_rad(angle, tolerance=5):
    '''
    Convert angle from degrees to radians with tolerance.

    **parameters:**
        angle: float
            Angle in degrees.
        tolerance: float
            Tolerance in degrees, default is 5.

    **returns:**
        bond_angle_rad: float
        Adjusted angle in radians.
    '''
    angle_variation = np.random.uniform(-tolerance, tolerance)
    adjusted_bond_angle = angle + angle_variation
    bond_angle_rad = np.deg2rad(adjusted_bond_angle)
    return bond_angle_rad


def find_key_or_value(key_or_value, dictionary):
    """
    Search for a key or value in a dictionary. If the key or value is
    found, returns the corresponding value or key.

    **parameters:**
        key_or_value : int or str
        dictionary : dict

    **returns:**
        The corresponding value or key if found, otherwise None.
    """
    if key_or_value in dictionary:
        return dictionary[key_or_value]

    for key, value in dictionary.items():
        if value == key_or_value:
            return key

    return None
# def add_dummy_to_point_of_extension(ase_atom, point_of_extension, mapped_indices):
#     '''
#     Add a dummy atom to the point of extension, ensuring that it forms an angle
#     of approximately 120° with the central atom and its neighbors, and is always
#     added outward.

#     parameters
#     ----------
#     ase_atom: Atoms object
#         The original atom structure.
#     point_of_extension: list
#         List of indices where dummy atoms are to be added.
#     mapped_indices: dict
#         Dictionary mapping internal indices to atomic indices.

#     returns
#     -------
#     new_atom: Atoms object
#         Atoms object containing the original atoms plus the newly added dummy atoms.
#     '''
#     distance = 1.0  # Distance for adding new dummy atom
#     all_positions = []

#     # Compute neighbor information
#     graph, bond_matrix = compute_ase_neighbour(ase_atom)

#     for idx in point_of_extension:
#         atom_idx = mapped_indices[idx]
#         neighbours = graph[atom_idx].tolist()
#         x0, y0, z0 = ase_atom[atom_idx].position

#         if len(neighbours) == 1:
#             # Add the dummy atom directly along the x-axis if there's only one neighbor
#             x = x0 + distance
#             y = y0
#             z = z0
#             new_position = [x, y, z]
#             all_positions.append(new_position)

#         elif len(neighbours) == 2:
#             # Get the positions of the two neighboring atoms
#             pos1 = ase_atom[neighbours[0]].position
#             pos2 = ase_atom[neighbours[1]].position

#             # Calculate vectors from the central atom to the neighbors
#             vec1 = pos1 - np.array([x0, y0, z0])
#             vec2 = pos2 - np.array([x0, y0, z0])

#             # Normalize the vectors
#             vec1 /= np.linalg.norm(vec1)
#             vec2 /= np.linalg.norm(vec2)

#             # Compute the normal to the plane (cross product of vec1 and vec2)
#             normal = np.cross(vec1, vec2)
#             normal /= np.linalg.norm(normal)

#             # Calculate the bisector direction (outward)
#             bisector = (vec1 + vec2) / np.linalg.norm(vec1 + vec2)

#             # Ensure the bisector points outward
#             if np.dot(bisector, normal) < 0:
#                 bisector = -bisector

#             # Rotate the bisector by 120 degrees to get the new atom's position
#             bond_angle_rad = angle_tolerance_to_rad(120)  # Desired angle

#             # Use the rotation matrix for in-plane rotation around the normal
#             cos_angle = np.cos(bond_angle_rad)
#             sin_angle = np.sin(bond_angle_rad)
#             rotation_matrix = np.array([
#                 [cos_angle + normal[0]**2 * (1 - cos_angle), normal[0] * normal[1] * (1 - cos_angle) - normal[2] * sin_angle, normal[0] * normal[2] * (1 - cos_angle) + normal[1] * sin_angle],
#                 [normal[1] * normal[0] * (1 - cos_angle) + normal[2] * sin_angle, cos_angle + normal[1]**2 * (1 - cos_angle), normal[1] * normal[2] * (1 - cos_angle) - normal[0] * sin_angle],
#                 [normal[2] * normal[0] * (1 - cos_angle) - normal[1] * sin_angle, normal[2] * normal[1] * (1 - cos_angle) + normal[0] * sin_angle, cos_angle + normal[2]**2 * (1 - cos_angle)]
#             ])

#             new_direction = np.dot(rotation_matrix, bisector)
#             new_position = np.array([x0, y0, z0]) + distance * new_direction

#             all_positions.append(new_position)

#     # Create a new Atoms object for the dummy atoms
#     symbols = ['X' for _ in all_positions]  # Use 'X' as a placeholder symbol for dummy atoms
#     new_atom = Atoms(symbols=symbols, positions=all_positions)

#     return new_atom


# def calculate_coordinates_atom(atom_type, bond_angle, atom_position, distance=0.9, tolerance=5):
#     """
#     Function to calculate to compute the coordinates of atoms to add to a give atoms. This
#     will mostly be used to add hydrogen atoms or dummies.
#     parameters
#     ----------
#     atom_type: str
#     bond_angle: float
#     atom_position: ndarray
#     distance: float
#     tolerance: float
#     returns the coordinates of the atoms

#     """
#     angle_variation = np.random.uniform(-tolerance, tolerance)
#     adjusted_bond_angle = bond_angle + angle_variation
#     bond_angle_rad = np.deg2rad(adjusted_bond_angle)

#     x0, y0, z0 = atom_position

#     if atom_type == "sp3":
#         # Tetrahedral: Place the hydrogen along the x-axis at distance d
#         x = x0 + distance
#         y = y0
#         z = z0
#     elif atom_type == "sp2":
#         # Trigonal planar: Place the hydrogen in the xy-plane
#         x = x0 + distance * np.cos(bond_angle_rad / 2)
#         y = y0 + distance * np.sin(bond_angle_rad / 2)
#         z = z0
#     elif atom_type == "sp":
#         # Linear: Place the hydrogen along the x-axis at distance d
#         x = x0 + distance
#         y = y0
#         z = z0
#     elif atom_type == "trigonal_pyramidal":
#         # Trigonal pyramidal: Place the hydrogen in the xz-plane
#         x = x0 + distance * np.cos(bond_angle_rad / 2)
#         y = y0
#         z = z0 + distance * np.sin(bond_angle_rad / 2)
#     elif atom_type == "bent":
#         # Bent: Place the hydrogen in the xz-plane
#         x = x0 + distance * np.cos(bond_angle_rad / 2)
#         y = y0
#         z = z0 + distance * np.sin(bond_angle_rad / 2)
#     else:
#         raise ValueError("Unsupported atom type")

#     return [x, y, z]
