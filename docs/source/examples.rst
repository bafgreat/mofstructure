How to Guide
============

This section provides step-by-step instructions on how to use `mofstructure` both from the command line and as a Python library. Whether you're new to this module or just need a refresher, this guide will help you get started with ease.

Run on the Command Line
==========================

One of the most powerful features of `mofstructure` is its ability to perform complex operations directly from the command line. Below, we walk you through how to deconstruct metal-organic frameworks (MOFs) into their building units and how to create a database of MOFs from multiple files.

Building Units
----------------

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

      connected_components, atoms_indices_at_breaking_point, porphyrin_checker, all_regions = mofdeconstructor.secondary_building_units(ase_atom)

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
