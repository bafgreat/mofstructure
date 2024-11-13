.. _updates-0.1.4:

Updates Version 0.1.4
======================

The new update enables the computation of open metal sites in CIFs. To use this functionality, run the following on the command line:

.. code-block:: bash

   mofstructure_database ciffolder --oms

Here, `ciffolder` corresponds to the directory/folder containing the CIF files.

After computation, the metal information will be found in a JSON file called `metal_info.json`. This file is located in the output folder that defaults to `MOFDb` if no other folder is provided.

.. note::

   Computing open metal sites is computationally expensive, especially if you intend to run it on a folder with many CIF files. If you are not interested in computing the open metal sites, simply run the command without the `--oms` option:

   .. code-block:: bash

      mofstructure_database ciffolder

   This command will generate a `MOFDb` folder without the `metal_info.json` file, but the code will run much faster.

Also, note that the `--oms` option is only provided for the `mofstructure_database` command and is not available for the `mofstructure` command, which targets a single CIF file. If you have a single CIF file and wish to compute open metal sites, simply put the CIF file in a folder and run the `mofstructure_database` command on that folder:

.. code-block:: bash

   mofstructure_database ciffolder --oms

   .. _updates-0.1.5:

Updates Version 0.1.5
======================

The new update enables users to include a RAD file when computing porosity using PyZeo. This allows users to specify the type of radii to use. If omitted, the default PyZeo radii will be used, which are covalent radii obtained from the CSD.

Currently, this functionality can only be used when using `mofstructure` as a library. This can be done as follows:

.. code-block:: python

   from mofstructure.porosity import zeo_calculation
   from ase.io import read

   ase_atom = read(filename)

   pore_data = zeo_calculation(ase_atom, rad_file='rad_file_name.rad')

.. note::

   `filename` is any ASE-readable crystal structure file, ideally a CIF file. Moreover, `rad_file_name.rad` is a file containing the radii of each element present in the structure file. This should be formatted as follows:

   .. code-block:: text

      element radii

   For example, for an MgO system, your RAD file should look like this:

   .. code-block:: text

      Mg 0.66
      O 1.84

   Also, note that if the radii file does not have the `.rad` extension like `rad_file_name.rad`, the default radii will be used.


.. _updates-0.1.6:

Updates Version 0.1.6
======================

Added new command line tools to expedite calculations, especially when working on a large database.

Compute Only Deconstruction
----------------------------

If you wish to only compute the deconstruction of MOFs without having to compute their porosity and open metal sites, then simply run the following command:

.. code-block:: bash

   mofstructure_building_units  cif_folder

Compute Only Porosity
----------------------

If you wish to only compute the porosity using default values (i.e., probe radius = 1.86, number of GCMC cycles = 10,000, and default CSD atomic radii), then run the following command:

.. code-block:: bash

   mofstructure_porosity cif_folder

However, if you wish to use another probe radius (e.g., 1.5) and GCMC cycles of 20,000, alongside custom atomic radii in a file called `rad.rad`, run the following command:

.. code-block:: bash

   mofstructure_porosity cif_folder -pr 1.5 -ns 20000 -rf rad.rad

Compute Only Open Metal Sites
------------------------------

If you are only interested in computing the open metal sites, then run the following command:

.. code-block:: bash

   mofstructure_oms cif_folder


.. _updates-0.1.7:

Updates Version 0.1.7
======================

1. Implemented a robust CI/CD using Git Actions.
2. Included `add_dummy` key to add dummy atoms to points of extension. This is important to effectively control the breaking point. These dummy atoms can then be replaced with hydrogen to fully neutralize the system.

.. note::

   Please don't use `add_dummy` when deconstructing to ligands and clusters. The `add_dummy` argument should be used only for SBUs, e.g.,

   .. code-block:: python

      connected_components, atoms_indices_at_breaking_point, porpyrin_checker, all_regions = MOF_deconstructor.secondary_building_units(ase_atom)
      metal_sbus, organic_sbus, building_unit_regions = MOF_deconstructor.find_unique_building_units(
          connected_components,
          atoms_indices_at_breaking_point,
          ase_atom,
          porpyrin_checker,
          all_regions,
          cheminfo=True,
          add_dummy=True
          )

      metal_sbus[0].write('test1.xyz')

Updates Version 0.1.8
======================
1. The main update here is to enable mofstructure to run on Python versions
3.9 to 3.13. We have done the neccessary tests but let us know if you have
any conflicts or bugs and we will fix it.


Updates Version 0.1.8.1
======================
Made rdkit to be an optional dependency so that mofstructure should
be compatible with Python 3.12, since there are no recent versions
of rdkit that are compatible with Python 3.12. Hence if you wish to use
rdkit, you should install it separately.

  .. code-block:: pip install rdkit

Updates Version 0.1.8.2
======================
Fixed the python dependency to be compatible with any Python 3.9 and above. 
  .. code-block:: pip install rdkit
