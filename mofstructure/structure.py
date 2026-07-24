#!/usr/bin/python
'''
High level entry point to mofstructure.

MOFstructure wraps a single crystal structure and exposes the analyses that are
normally wanted together: guest removal, porosity, deconstruction into building
units, topology and open metal sites. Everything is computed on the guest free
system, so a structure only has to be read once.

    from mofstructure import structure

    mof = structure.MOFstructure(filename='UiO-66.cif')
    pores = mof.get_porosity()
    metal_sbus, organic_sbus = mof.get_sbu()

The individual algorithms live in mofdeconstructor, porosity and systre, and can
be called directly when more control is needed.
'''
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
import os
from functools import reduce
import operator
import argparse
import pandas as pd
import shutil
import tempfile
import json
import numpy as np
import logging
from ase.io import read
from omsdetector_forked import MofCollection, mof
import mofstructure.mofdeconstructor as MOF_deconstructor
from mofstructure.porosity import zeo_calculation
import mofstructure.filetyper as read_write
from mofstructure.systre import SystreTopology

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


class MOFstructure(object):
    '''
    A single framework and the analyses that can be run on it.

    Give it either an ASE atoms object or a path to any file ASE can read. The
    methods return their results rather than writing them, so they can be mixed
    freely:

        mof = MOFstructure(filename='UiO-66.cif')
        pores = mof.get_porosity()
        metal_sbus, organic_sbus = mof.get_sbu()
        topology = mof.get_topology()

    Guest molecules are stripped internally before each analysis, so a structure
    that still contains solvent does not need cleaning up first.

    **parameters:**
        ase_atoms: ASE atoms object, used in preference to filename
        filename: path to a cif or any other ASE readable structure file
    '''

    def __init__(self,
                 ase_atoms=None,
                 filename=None
                 ):
        """
        A global method for extracting useful information
        """
        if ase_atoms is not None:
            self.ase_atoms = ase_atoms
        else:
            self.ase_atoms = read(filename)

    def remove_guest(self):
        '''
        Simple function to remove guest molecules in porous system.
        Note that this can work for any periodic system.

        **return:**
            ase_atoms (ase.Atoms): atom object with guest removed.
        '''
        index_non_guest = MOF_deconstructor.remove_unbound_guest(self.ase_atoms)
        return self.ase_atoms[index_non_guest]

    def get_sbu(self, wrap_system=True, cheminfo=True, add_dummy=False):
        """
        Extract the metal and linker secondary building units from the guest free system.

        **parameter**
            wrap_system (bool): If True, removes the effects of periodicity by merging the system.
            cheminfo (bool): If True, computes cheminformatic identifiers such as SMILES, InChI, and InChIKey.
            add_dummy (bool): If True, adds dummy atoms at the points of extension.

        **return:**
            metal_sbu (list): A list of unique metal secondary building units
            linker_sbu (list): A list of unique organic secondary building units
        """
        guest_free_atoms = self.remove_guest()
        connected_components, atoms_indices_at_breaking_point, porpyrin_checker, all_regions, _ = MOF_deconstructor.secondary_building_units(guest_free_atoms)
        if len(connected_components) > 0:
            metal_sbus, organic_sbus, _ = MOF_deconstructor.\
                find_unique_building_units(connected_components,
                                           atoms_indices_at_breaking_point,
                                           guest_free_atoms,
                                           porpyrin_checker,
                                           all_regions,
                                           wrap_system=wrap_system,
                                           cheminfo=cheminfo,
                                           add_dummy=add_dummy
                                           )
            return metal_sbus, organic_sbus
        else:
            logging.error(
                "Sorry, we were unable to successfully deconstruct your system."
                "Send us an email so that we can implement better rules for"
                "deconstructing your system."
            )
            return None

    def get_ligands(self, wrap_system=True, cheminfo=True, add_dummy=False):
        """
        Extract the metal and linker secondary building units from the guest free system.

        **parameter:**
            wrap_system (bool): If True, removes the effects of periodicity by merging the system.
            cheminfo (bool): If True, computes cheminformatic identifiers such as SMILES, InChI, and InChIKey.
            add_dummy (bool): If True, adds dummy atoms at the points of extension.

        **return:**
            metal_clusters (list): A list of unique metal atoms or clusters.
            organic_ligands (list): A list of unique organic ligands.
        """
        guest_free_atoms = self.remove_guest()
        connected_components, atoms_indices_at_breaking_point, porpyrin_checker, all_regions, _ = MOF_deconstructor.ligands_and_metal_clusters(guest_free_atoms)
        if len(connected_components) > 0:
            metal_clusters, organic_ligands, _ = MOF_deconstructor.\
                find_unique_building_units(connected_components,
                                           atoms_indices_at_breaking_point,
                                           guest_free_atoms,
                                           porpyrin_checker,
                                           all_regions,
                                           wrap_system=wrap_system,
                                           cheminfo=cheminfo,
                                           add_dummy=add_dummy
                                           )
            return metal_clusters, organic_ligands
        else:
            logging.error(
                "Sorry, we were unable to successfully deconstruct your system."
                "Send us an email so that we can implement better rules for"
                "deconstructing your system."
            )
            return None

    def get_porosity(self,
                     probe_radius=1.86,
                     number_of_steps=10000,
                     rad_file=None,
                     high_accuracy=True):
        '''
        A function to compute porosity data for a system.

        **parameters**

            probe_radius (float): Radius of the probe (default: 1.86).
            number_of_steps (int): Number of GCMC simulation cycles (default: 10000).
            high_accuracy (bool): If True, perform high-accuracy computations.
            rad_file: Optional file containing user defined atom radii. Must have the `.rad` extension

        **return:**
            pore (dict): A dictionary containing:
                - AV_Volume_fraction: Accessible volume void fraction.
                - AV_A^3: Accessible volume in A^2.
                - AV_cm^3/g: Accessible volume in cm³/g. This value is often infinite because it divides the computed
                    volume by Avogadro's number.
                - ASA_A^2: Accessible surface area in A^2.
                - ASA_m^2/cm^3: Accessible surface area in m2/cm3.
                - Number_of_channels: Number of channels (i.e., pores) present in the system.
                - LCD_A: The largest cavity diameter, defined as the diameter of the largest sphere that can be
                    inserted into the porous system without overlapping any atoms.
                - lfpd_A: The largest included sphere along the free sphere path, i.e., the largest sphere that can be
                    inserted into the pore.
                - PLD_A: The pore limiting diameter, defined as the largest sphere that can freely diffuse through the
                    porous network without overlapping any atoms.
        '''
        guest_free_atoms = self.remove_guest()
        if rad_file is None:
            pores = zeo_calculation(guest_free_atoms,
                                    probe_radius=probe_radius,
                                    number_of_steps=number_of_steps,
                                    high_accuracy=high_accuracy
                                    )
        else:
            pores = zeo_calculation(guest_free_atoms,
                                    probe_radius=probe_radius,
                                    number_of_steps=number_of_steps,
                                    high_accuracy=high_accuracy,
                                    rad_file=rad_file)
        return read_write.convert_numpy_types(pores)

    def get_topology(
        self,
        method="all_node",
        *,
        decimals=8,
        include_edge_centers=True,
        fallback_to_input_cgd=False,
    ):
        """
        Compute topology information for the guest-free system.

        **parameters:**
            method: str
                Topology extraction method passed to SystreTopology.

            decimals: int
                Number of decimal places used when hashing the relaxed topology payload.

            include_edge_centers: bool
                If True, include edge-center comments in the CRYSTAL2 text.

            fallback_to_input_cgd: bool
                If True, return a CRYSTAL2 wrapper from the input CGD when no
                relaxed component can be parsed from Systre output.

        **return:**
            python dictionary
                Mapping containing:
                    - topology
                    - dimension
                    - td10
                    - topology_hash
                    - cgd
        """
        guest_free_atoms = self.remove_guest()

        runner = SystreTopology(
            guest_free_atoms,
            method=method,
            name="net",
            keep_tmp=False,
        )

        result = runner.identify()
        comp = runner.best_component()

        if comp is None:
            return {
                "topology": result.topology,
                "dimension": None,
                "td10": None,
                "topology_hash": None,
                "cgd": runner.crystal2_text(
                    include_edge_centers=include_edge_centers,
                    fallback_to_input_cgd=fallback_to_input_cgd,
                ),
            }

        return {
            "topology": result.topology,
            "dimension": comp.dimension,
            "td10": comp.td10,
            "topology_hash": comp.topology_hash(decimals=decimals),
            "cgd": comp.crystal2_text(
                name=comp.rcsr_name or "net",
                include_edge_centers=include_edge_centers,
            ),
        }

    def get_oms(self):
        """
        Function to compute open metal sites
        """
        general_info = {}
        overlap = MOF_deconstructor.inter_atomic_distance_check(self.ase_atoms)
        general_info['has_overlapping_atoms'] = not overlap

        if len(self.ase_atoms) > 3000:
            print(f'The system size is very large: {len(self.ase_atoms)} atoms')
            print('Hence oms computation will take a while ')
            print('Thanks for your patience!!!')

        lattice = self.ase_atoms.get_cell().tolist()
        species = self.ase_atoms.get_chemical_symbols()
        coords = self.ase_atoms.get_positions()
        oms_structure = mof.MofStructure(lattice=lattice,
                                         species=species,
                                         coords=coords,
                                         coords_are_cartesian=True,
                                         name='oms_system'
                                         )

        coordination_sphere = oms_structure.metal_coord_spheres
        tmp_dir = tempfile.mkdtemp()
        try:
            oms_structure.analyze_metals(tmp_dir)

            json_files = [f for f in os.listdir(tmp_dir) if f.endswith('.json')]
            if not json_files:
                raise FileNotFoundError(f"No JSON file found in temporary directory: {tmp_dir}")


            json_file = os.path.join(tmp_dir, json_files[0])
            data = read_write.load_data(json_file)


            general_info["metals"] = data.get("metal_species")
            general_info["has_oms"] = data.get("has_oms")
            general_info["density"] = data.get("density")
            general_info["uc_volume"] = data.get("uc_volume")
            general_info["error_in_systems"] = data.get("problematic")


            metal_sites_dic = data.get("metal_sites")
            if metal_sites_dic is not None:
                unique_tuples = {
                    (d["metal"], d["number_of_linkers"], d["is_open"])
                    for d in metal_sites_dic
                }
                unique_metal = [
                    {
                        "metal": metal,
                        "coordination_number": number,
                        "is_open": is_open,
                        "enviroment": get_coordination_environment(coordination_sphere, metal, number)
                    }
                    for metal, number, is_open in unique_tuples
                ]
                oms_metal = [i.get("metal") for i in unique_metal if i["is_open"]]
                general_info["metal_info"] = unique_metal
                general_info["open_metals"] = oms_metal
        finally:
            shutil.rmtree(tmp_dir)


        return general_info



def get_coordination_environment(coordination_spheres, metal_element, coordination_number):
    """
    Returns the coordination environment (neighbors) of a specified metal element
    only if the number of neighbors equals the given coordination number.

    **parameters:**
        coordination_spheres (list): A list of metal coordination spheres (MetalSite objects).
        metal_element (str): The metal element symbol (e.g., 'Ag').
        coordination_number (int): The desired number of neighbors (coordination number).

    Returns:
        List[str]: A list of species symbols (neighbors) bonded to the metal center.
                   Returns an empty list if no matching environment is found.
    """
    metal_element = str(metal_element)

    for metal_site in coordination_spheres:
        # Convert the species list to strings.
        species_str = [str(s) for s in metal_site.species]

        # Early check: if the total number of atoms is not metal + neighbors, skip.
        if len(species_str) != coordination_number + 1:
            continue

        try:
            center_idx = species_str.index(metal_element)
        except ValueError:
            continue

        # Get the neighbors by removing the metal element.
        neighbors = species_str[:center_idx] + species_str[center_idx+1:]
        if len(neighbors) == coordination_number:
            return neighbors

    return []
