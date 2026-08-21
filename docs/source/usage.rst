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

The database will be saved in a folder named `MOFstructureDB` in the current directory.

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

The record always holds the same keys - ``av_volume_fraction``, ``av_a3``,
``asa_a2``, ``asa_m2_per_cm3``, ``number_of_channels``, ``lcd_a``, ``pld_a``,
``lfpd_a`` and ``porosity_status``. A structure zeo++ cannot handle, or one
that runs past ``timeout`` seconds, comes back with ``None`` in place of each
number and the reason in ``porosity_status``, so a batch never stalls and a
table never acquires holes.

Example 4: Deconstructing MOFs into SBUs and linkers
-----------------------------------------------------

To identify and extract SBUs and linkers from the MOF:

.. code-block:: python

   connected_components, atoms_indices_at_breaking_point, porpyrin_checker, all_regions, breaking_pairs = MOF_deconstructor.secondary_building_units(ase_atom)

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
To determine the topology of a framework. The kind of framework is worked out
from the structure, so MOFs, COFs and zeolites can be passed alike:

.. code-block:: python

   from ase.io import read
   from mofstructure.topology import analyse, analyse_methods

   # 1) From a CGD periodic graph
   print(analyse("net.cgd")["topology"])

   # 2) From a CIF, naming the deconstruction
   print(analyse("UiO-66.cif", method="all_node")["topology"])

   # 3) From ASE Atoms, letting it choose
   print(analyse(read("UiO-66.cif"))["topology"])

   # 4) A MOF has more than one defensible net; report them all
   for method, record in analyse_methods("UiO-66.cif").items():
       print(method, record["topology"])

The record also carries ``key``, the canonical key. It is unique to the net and
unchanged by supercell, atom order or origin, so two structures with the same
key have the same topology whether or not any archive names it.



Example 7: Computing topology from the command line
-----------------------------------------------------
To compute the topology of a MOF from the command line, run:

.. code-block:: bash

   mofstructure_topology cif_folder

This computes the topology of every CIF file in the folder and appends the
records to ``MOFstructureDB/Structure_Data/topology_data.json`` - the same file
``mofstructure_database`` writes to, so one folder holds one structure
database however it was built. A ``topology_data.csv`` summary is written
beside it, one row per structure. Use ``-s`` to name a different directory,
``--json`` to write the full records to a file of your choosing, and
``--no-save`` to print without writing.

A single file, several files or a folder may be named. Each structure is
printed as it finishes and a tally closes the run; ``--quiet`` keeps the tally
only, and ``-v`` reports each deconstruction as it is built.

A MOF admits more than one node definition. Pass ``--method`` to pick one of
``sbus``, ``all_node``, ``single_node`` or ``ligand_cluster``, or
``--all-methods`` to record each of them under ``<structure>:<method>``.


Example 8: Writing the net as a CGD file
-----------------------------------------------------
``mofstructure_generate_cgd`` writes the net itself as a CGD
``PERIODIC_GRAPH``. ``--method`` takes the same node definitions, and ``auto``
reads the material from the structure:

.. code-block:: bash

   mofstructure_generate_cgd HKUST-1.cif
   mofstructure_generate_cgd zeolite.cif --method zeol -o zeolite.cgd

``--embedding`` chooses the geometry written for the net.
``deconstruction``, the default, writes the crystal's own coordinates and needs
no identification. ``ideal`` writes the canonical barycentric placement, with
the cell fixed by the symmetry of the net alone; it is reproducible to the
digit and carries the canonical key, so an unnamed net stays identifiable from
its own file. ``refined`` relaxes that towards edges of equal length, which is
what a tool building on the net needs.

.. code-block:: bash

   mofstructure_generate_cgd HKUST-1.cif --embedding ideal
   mofstructure_generate_cgd UiO-66.cif --embedding refined

The barycentric placement minimises the sum of squared edge lengths, so real
nets emerge with the longest edge two or three times the shortest. Refining
typically brings that spread from 2.2 down to 1.0, but its result depends on
the optimiser and is not reproducible to the last digit, so every file records
which embedding it holds. Where refining cannot beat the exact placement, the
exact one is written and the command says so.