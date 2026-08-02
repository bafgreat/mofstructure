How to Guide
============
This section provides detailed
step-by-step examples demonstrating how to use `mofstructure`
both as a Python library and from the command line.
The aim is to guide users from a raw crystal structure
(e.g. CIF file) to meaningful structural, topological and chemical insights.


Quick Start Guide
===================

The `MOFstructure` class provides a simple
and unified interface to most functionalities

.. code-block:: python

   from mofstructure import structure
   import mofstructure.filetyper as read_write

   # Read a CIF file
   cif_file = 'path_to_your_cif_file.cif'
   mofdata = structure.MOFStructure(filename=cif_file)
   # can also directly parse an ase_atoms object
   # mofdata = structure.MOFStructure(ase_atoms=ase_atoms)

   # remove unbound guest molecules
   guest_free_mof = mofdata.remove_guest()

   # get the topology of the MOF
   topology = mofdata.get_topology(method="all_node") # "all_node", "single_node" or "sbus"

   print(topology['topology']) # prints the topology name, e.g. pcu, dia, etc.
   print(topology['dimension']) # prints the dimension of the topology, e.g. 3 for 3D, 2 for 2D, etc.
   print(topology['td10']) # prints the td10 value of the topology, which is a measure of the
   # get metal and organic SBUs
   # This will return a list of ase_atoms objects
   metal_sbus, organic_sbus = mofdata.get_sbu()

   # get organic ligands
   _, organic_ligands = mofdata.get_ligands()

   # porosity
   pores = mofdata.get_porosity(probe_radius=1.86, number_of_steps=10000, rad_file=None, high_accuracy=True)

   # get open metal sites
   open_metal_sites = mofdata.get_oms()

Reading Structures
====================
`mofstructure` supports all formats readable by ASE (e.g. CIF, POSCAR, XYZ).

.. code-block:: python

   from mofstructure import structure

   mof = structure.MOFstructure(filename="structure.cif")

Alternatively an ASE Atoms object can be passed directly:

.. code-block:: python

   from mofstructure import structure
   mof = structure.MOFstructure(ase_atoms=ase_atoms)


Removing Guest Molecules
========================
Many experimentally resolved MOFs contain solvent
or guest molecules that are not part of the framework.
These must be removed before analysis.

.. code-block:: python

   clean_structure = mof.remove_guest()

This codes identifies disconnected components in the periodic graph and removes unbound fragments.

Deconstruction: Building Units
===============================
A central feature of mofstructure is the decomposition of a framework into:
- Metal secondary building units (SBUs)
- Organic SBUs (linkers)
- Organic ligands
- metal clusters

.. code-block:: python

   metal_sbus, organic_sbus = mof.get_sbu()
   metal_clusters, organic_ligands = mof.get_ligands()

Each building unit is returned as an ASE Atoms object with additional metadata stored in .info.

.. code-block:: python

   for sbu in metal_sbus:
      print(sbu.info["sbu_type"])
      print(sbu.info["inchikey"])

Available information includes:
- SMILES
- InChI / InChIKey
- SBU type (e.g. paddlewheel, rod-like)
- Points of extension

Saving building units:
-----------------------

.. code-block:: python

   for i, sbu in enumerate(metal_sbus):
      sbu.write(f"metal_sbu_{i}.cif")

Topology Analysis
=================

Topology is determined using a graph representation of the framework and analysed using Systre.
This is possible becuase of the deconstruction into building units, which allows for different levels of abstraction.

.. code-block:: python

   topology = mof.get_topology(method="all_node")
   print(topology["topology"]) # e.g. pcu, dia
   print(topology["dimension"]) # 2 or 3
   print(topology["td10"])

Available methods:

- all_node: rod SBUs split into their atoms (CrystalNets AllNodes)
- single_node: all_node with organic groups merged (CrystalNets SingleNodes)
- sbus: each SBU as one node

The output includes:
- RCSR topology name
- Dimensionality
- Topological descriptors
- TD10 value
- Topology hash for uniqueness
- systre optimised cgd string representation

Porosity Analysis
=================
Porosity properties are computed using our python wrapper around Zeo++ called pyzeo.

.. code-block:: python

   pores = mof.get_porosity(
   probe_radius=1.86,
   number_of_steps=10000,
   high_accuracy=True
   )
   print(pores["AV_Volume_fraction"])
   print(pores["ASA_m2_cm3"])

Typical outputs include:
- Accessible volume fraction
- Accessible surface area
- Pore limiting diameter (PLD)
- Largest cavity diameter (LCD)

For high-quality results, `high_accuracy=True` is recommended.


Open Metal Sites (OMS)
=======================
Open metal sites are important for adsorption and catalysis.

   .. code-block:: python

   oms = mof.get_oms()

The returned data typically includes:
- Metal identity
- Coordination environment
- Atomic indices

Run on the Command Line
==========================

One of the most powerful features of `mofstructure` is its ability to perform complex operations directly from the command line. Below, we walk you through how to deconstruct metal-organic frameworks (MOFs) into their building units and how to create a database of MOFs from multiple files.

Deconstruct a structure:
------------------------

If you have a CIF file (or any file format that ASE can read, such as POSCAR, XYZ, etc.) containing a MOF, you can deconstruct it into its constituent building units using a simple command. This command processes the MOF structure and saves the results in an organized folder structure.

To begin, navigate to the directory containing your CIF file and execute the following command:

.. code-block:: bash

   mofstructure cif_file

Here, `cif_file` should be replaced with the name of your actual CIF file. The script will automatically deconstruct the MOF present in the file and generate the output in a folder named `MOF_building_units` within the current directory.

Custom Output Directory
------------------------

If you wish to store the results in a specific directory, rather than the default `MOF_building_units` folder, you can specify the path to your desired output folder like this:

.. code-block:: bash

   mofstructure cif_file path_to_result_folder

Replace `path_to_result_folder` with the full or relative path to the directory where you want the output saved.

Processing Multiple CIF Files
==============================

For cases where you have multiple CIF files that need to be processed, you can automate the process by running a loop in Python or a shell script. The results for each file will be saved in the specified directory:

.. code-block:: python

   for cif in ciffiles:
       mofstructure cif path_to_result

In this example, `ciffiles` is a list of all CIF file paths that you want to process. The script will iterate over each file, deconstruct the MOF, and save the output accordingly.

Creating a Database
---------------------

If you have a collection of CIF files stored in a single directory, you can easily create a comprehensive database of MOFs. This database will compile all the MOF structures into a neatly organized format, making it easier to manage and analyze large datasets.

To create the database, run the following command:

.. code-block:: bash

   mofstructure_database ciffolder

Here, `ciffolder` should be replaced with the path to the folder containing all your CIF files. The output will be automatically saved in a folder named `MOFDb` within your current working directory.

Custom Database Output Directory
---------------------------------

If you prefer to save the database in a different location, you can specify the desired output path directly in the command:

.. code-block:: bash

   mofstructure_database ciffolder path_to_result

Replace `path_to_result` with the path to the folder where you want the database to be stored. This flexibility allows you to organize your work according to your preferences.

Use as a Library
================

In addition to command-line usage, `mofstructure` can also be used as a Python library, providing more granular control and flexibility for advanced users. Below are the key steps to get started with `mofstructure` as a library.

1. **Importing the Module**

   Begin by importing the necessary components from `mofstructure` and any other dependencies required for your workflow:

   .. code-block:: python

      from mofstructure import mofdeconstructor
      from mofstructure import porosity
      from mofstructure import buildingunits
      from ase.io import read, write
      import pandas as pd

2. **Reading a MOF File Using ASE**

   Use the ASE (Atomic Simulation Environment) library to read the CIF file or any other supported file format:

   .. code-block:: python

      ase_atom = read(cif_file)

   Here, `cif_file` is the path to your CIF file. The `read` function loads the structure into an `ase.Atoms` object, which can then be manipulated using `mofstructure`.

3. **Removal of Unbound Guest Molecules**

   If your MOF structure contains unbound guest molecules, you can easily remove them using the following command:

   .. image:: images/guest_removal.gif
      :alt: Guest Removal

   .. code-block:: python

      no_guest_indices = mofdeconstructor.remove_unbound_guest(ase_atom)
      no_guest_atom = ase_atom[no_guest_indices]

   The `remove_unbound_guest` function returns the indices of atoms that are not part of unbound guest molecules, allowing you to filter them out and work with a cleaner structure.

4. **Computing Porosity**

   To compute porosity properties such as pore size distribution, you can use the `zeo_calculation` function:

   .. code-block:: python

      pores = porosity.zeo_calculation(ase_atom, probe_radius=1.86, number_of_steps=5000)
      df = pd.DataFrame(pores, index=[0])
      df.to_csv('pore.csv')

   This command performs a porosity analysis using a probe with a specified radius and saves the results in a CSV file.

5. **Identifying SBUs and Ligands**

   Deconstruct the MOF into its Secondary Building Units (SBUs) and ligands:

   .. image:: images/deconstruction.gif
      :alt: MOF Deconstruction

   .. code-block:: python

      connected_components, atoms_indices_at_breaking_point, porphyrin_checker, all_regions, breaking_pairs = mofdeconstructor.secondary_building_units(ase_atom)

      metal_sbus, organic_sbus, building_unit_regions = mofdeconstructor.find_unique_building_units(
          connected_components,
          atoms_indices_at_breaking_point,
          ase_atom,
          porphyrin_checker,
          all_regions,
          cheminfo=True

      )

   By setting `cheminfo=True`, `mofstructure` calls Open Babel to compute cheminformatic information such as SMILES, InChI, and InChIKey, which are stored in `ase_atom.info`. The `metal_sbus` and `organic_sbus` lists contain all unique instances of the metal and organic SBUs, respectively.

Extracting Cheminformatic Information
--------------------------------------

To access and save the cheminformatic data for each SBU, you can iterate through the list of building units as shown below:

.. code-block:: python

   for i, sbu in enumerate(metal_sbus):
       smi = sbu.info['smi']
       inchi = sbu.info['inchi']
       inchikey = sbu.info['inchikey']
       # For SBUs only
       number_of_points_of_extension = sbu.info['points_of_extension']
       # For Metal SBUs only
       sbu_type = sbu.info['sbu_type']  # sbu_type could be rodlike, IRMOF, UIO66, paddlewheel, etc.
       # Save the SBU structure to a file
       sbu.write('metal_sbu_'+str(i)+'.cif')

This code snippet extracts relevant cheminformatic information for each SBU and saves the SBU structures in separate CIF files, named sequentially according to their index.

By following this comprehensive guide, you should now be well-equipped to start using `mofstructure` for your MOF research and projects. Whether you are manipulating structures from the command line or within a Python script, `mofstructure` offers the flexibility and power to handle a wide range of tasks in MOF analysis.
