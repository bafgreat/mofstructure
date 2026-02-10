Quick Start Guide
=================

This Quick Start Guide provides simple examples to help you begin using `mofstructure` for working with metal-organic frameworks (MOFs). We'll cover basic command-line usage and how to use `mofstructure` as a Python library.

Running on the Command Line
==================================

Example 1: Deconstructing a MOF into Building Units
-------------------------------------------------------

To deconstruct a MOF from a CIF file into its building units, run the following command:

.. code-block:: bash

   mofstructure example_mof.cif

This command processes the `example_mof.cif` file and saves the output in a folder named `MOF_building_units` within the current directory.

Example 2: Specifying an Output Directory
-----------------------------------------------

If you want to save the deconstructed MOF in a specific folder, use this command:

.. code-block:: bash

   mofstructure example_mof.cif /path/to/output_folder

Replace `/path/to/output_folder` with your desired output directory.

Example 3: Creating a Database from Multiple CIF Files
---------------------------------------------------------

To create a database from a folder containing multiple CIF files, use the following command:

.. code-block:: bash

   mofstructure_database /path/to/cif_folder

The database will be saved in a folder named `MOFDb` in the current directory.

Using `mofstructure` as a Library
==========================================

Example 1: Importing the Module and Reading a CIF File
-------------------------------------------------------

Start by importing the necessary modules and reading a CIF file using ASE:

.. code-block:: python

   from mofstructure import mofdeconstructor
   from ase.io import read

   ase_atom = read('example_mof.cif')

Example 2: Removing Unbound Guest Molecules
________________________________________________

To remove unbound guest molecules from the structure, use:

.. code-block:: python

   no_guest_indices = mofdeconstructor.remove_unbound_guest(ase_atom)
   no_guest_atom = ase_atom[no_guest_indices]

Example 3: Computing Porosity
------------------------------

To compute the porosity of the MOF, run:

.. code-block:: python

   from mofstructure import porosity

   pores = porosity.zeo_calculation(ase_atom, probe_radius=1.86)
   print(pores)

Example 4: Deconstructing MOFs into SBUs and linkers
-----------------------------------------------------

To identify and extract SBUs and linkers from the MOF:

.. code-block:: python

   connected_components, atoms_indices_at_breaking_point, porpyrin_checker, all_regions = MOF_deconstructor.secondary_building_units(ase_atom)

    metal_sbus, organic_sbus, _ = MOF_deconstructor.find_unique_building_units(
                connected_components,
                atoms_indices_at_breaking_point,
                ase_atom, porpyrin_checker,
                all_regions,
                cheminfo=True,
                add_dummy=True
            )

This code will output the SBUs and linkers along with their cheminformatic information.


Example 6: Determining the Topology of a MOF
-----------------------------------------------------
To determine the topology of a MOF using Systre:
.. code-block:: python

   from mofstructure.systre import identify_topology

   # 1) From a CGD file
   res = identify_topology("net.cgd", input_is_cgd=True)
   print(res.topology)

   # 2) From a CIF (generate CGD then run systre)
   res = identify_topology("UiO-66.cif", method="all_node")
   print(res.topology)

   # 3) From ASE Atoms
   atoms = read("UiO-66.cif")
   res = identify_topology(atoms, method="sbus")
   print(res.topology)



Example 7: Computing topoology from the command line
-----------------------------------------------------
To compute the topology of a MOF from the command line, run:
.. code-block:: bash

   mofstructure_topology cif_folder
This command will compute the topology of all CIF files in the specified folder and save the results in a file named `topology_results.csv` in the current directory.