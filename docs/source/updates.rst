Updates Version 0.1.8.9
=======================

This release corrects two topology bugs, fixes a crash when adding dummy atoms,
and removes a redundant option. Deconstruction was already correct in the
topology cases; only the net handed to Systre was wrong.

Topology
--------

**Polytopic linkers were contracted into a clique.** Every linker was collapsed
into edges directly between the metals it touched. For a ditopic linker that is
one edge, which is right, but a tritopic linker (BTC and higher) became a
triangle of edges, inflating each metal's connectivity. HKUST-1, a known
``tbo`` net, came out ``reo``. A polytopic linker is now its own branch-point
node joined to each metal it bridges, so HKUST-1 gives ``tbo``. Ditopic
frameworks are unaffected: UiO-66 stays ``fcu``, MOF-5 and pillared paddlewheels
stay ``pcu``.

**Rod SBUs lost their chain connectivity.** An infinite metal-oxo rod (MIL-53
and similar) is periodic within itself, but that periodicity was discarded when
the rod was contracted to a node, collapsing the net to a two-dimensional
``sql``. ``method="all_node"`` now splits a rod into its atoms -- each metal and
each bridging carboxyl carbon becomes a node -- recovering the true net. MIL-53
gives ``rna`` (matching CrystalNets AllNodes and mofid), while ``method="sbus"``
keeps the rod as one node (``pcu``). Discrete SBUs are unchanged (HKUST-1
``tbo``, UiO-66 ``fcu``).

**New ``method="single_node"``.** This is the coarsening CrystalNets calls
SingleNodes: the all-node net with each connected group of organic vertices
merged to one vertex, metals left separate. MIL-53 gives ``bpq``, matching
CrystalNets. ``mofstructure_topology --method all`` now records sbus, all_node
and single_node together.

Fixes
-----

``find_unique_building_units(..., add_dummy=True)`` crashed with an
``IndexError`` because the helper that locates the far side of a broken bond
only matched the first atom of the pair. Matching is now symmetric, so a dummy
atom is placed at every point of extension.

Command line
------------

``mofstructure_topology`` gained ``--method all``, which runs the node
definitions (``sbus``, ``all_node``, ``single_node``) and records the nets in
one entry per structure. The JSON nests them under a ``topologies`` key and the
CSV gives each its own columns, ready to load into a database.

The ``ligand_cluster`` method was removed: it duplicated ``sbus`` on simple
frameworks and gave non-standard nets on the rest. Use ``sbus``, ``all_node``
or ``single_node``.

The console scripts were made consistent: every script accepts
``-v/--verbose``, ``mofstructure_database`` accepts ``--method`` (with
``--topology_method`` kept as an alias), and ``-o`` now means output
everywhere. ``mofstructure_database`` previously used ``-o`` for ``--oms``;
that short form was removed, so use ``--oms``.

Changes
-------

The ``connect_mode`` argument was removed from ``build_cgd`` and the
``mofstructure_topology`` CLI. Linker contraction no longer has a clique/chain
choice: ditopic linkers become edges and polytopic linkers become nodes.

Updates Version 0.1.8.7
=======================

This release fixes ligand naming, which never resolved before, and stops a
crash in the porosity code that could terminate a database run.

IUPAC ligand names now resolve
------------------------------

``ligand_names`` previously came back as ``null`` for every structure. Two
separate problems had to be fixed before a name could match.

The shipped name database was scraped from PubChem, whose SMILES come from the
CACTVS toolkit. Canonical SMILES are only canonical within one toolkit, so a
CACTVS string never equals the OpenBabel string that ``compute_smi`` produces.

A ligand cut out of a framework also carries dangling valences at its points of
extension. Terephthalate leaves deconstruction two hydrogens short of the
terephthalic acid PubChem stores, so the formula, the SMILES and even the
InChIKey connectivity block all differ.

Both sides are now normalised the same way. Open valences are saturated to
recover the neutral parent molecule, and every molecule is indexed under its
full InChIKey, its canonical SMILES and its InChIKey connectivity block.

``get_ligands()`` returns ASE atoms objects, not names. Look the name up from
the SMILES carried on each fragment:

.. code-block:: python

   from mofstructure import structure
   from mofstructure.filetyper import load_iupac_names
   from mofstructure.mofdeconstructor import lookup_iupac_name

   iupac_names = load_iupac_names()

   mof = structure.MOFstructure(filename='RUBTAK01.cif')
   _, ligands = mof.get_ligands()

   for ligand in ligands:
       print(lookup_iupac_name(ligand.info['smi'], iupac_names))
   # terephthalic acid

The command line tools do this for you and write the result to the
``ligand_names`` field of ``ligands_data.json``, which previously held only
``null``.

Porosity no longer terminates a run
------------------------------------

zeo++ calls ``abort()`` when its Voronoi decomposition fails an internal volume
check, which raises SIGABRT and cannot be caught by any ``try``/``except``. A
single awkward framework previously killed the whole process.

``zeo_calculation`` now runs in a child interpreter and returns an empty
dictionary when the child dies, so such structures are recorded as having no
porosity data instead of ending the job.

Topology in the database workflow
----------------------------------

``mofstructure_database`` accepts ``-t/--topology``, writing
``topology_data.json`` and a CSV summary. It is off by default because Systre
runs on the JVM and costs a few seconds per structure.

.. code-block:: bash

   mofstructure_database cif_folder -t
   mofstructure_database cif_folder -t --topology_method all_node

Other fixes
-----------

- The CSV summaries no longer fail when a structure was recorded as ``null``.
- ``mofstructure_generate_cgd`` pointed at a module that does not exist and
  failed on every install. It now resolves correctly.
- Module level documentation across the package, and a docstring for the
  ``MOFstructure`` class.

Updates Version 0.1.8.6
=======================

Major update in the code structuring with
 the inclusion of TopologyExtractor class which makes
 the entire code more modular and easy to use.
 In addition, the code can now guess the names of common ligands.
 Moreover, we have added a new command line tool to compute the topology of MOFs using Systre.
 This tool is called `mofstructure_topology` and can be used as follows:

.. code-block:: bash

   mofstructure_topology cif_folder

   mofstructure_topology cif_folder --method all_node
   mofstructure_topology cif_folder --method sbus

   # Can be directly used on a single CIF file as well
   mofstructure_topology cif_file.cif


The topology can also be computing as a library as follows:

.. code-block:: python

   from mofstructure.systre import identify_topology
   from ase.io import read

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


It can also be called from structure as follows:

.. code-block:: python

   from mofstructure import structure
   from ase.io import read

   atoms = read("UiO-66.cif")
   mofdata = structure.MOFStructure(ase_atoms=atoms)
   topology = mofdata.get_topology(method="sbus")
   print(topology)

Updates Version 0.1.8.3
=======================
Major update in the code structuring with the inclusion of MOFstructure class
which makes the entire code more modular and easy to use. In addition,
the code can now guess the names of common ligands.

Updates Version 0.1.8.2
=======================
Fixed the python dependency to be compatible with any Python 3.9 and above.

  .. code-block:: bash

     pip install rdkit

Updates Version 0.1.8.1
=======================
Made rdkit to be an optional dependency so that mofstructure should
be compatible with Python 3.12, since there are no recent versions
of rdkit that are compatible with Python 3.12. Hence if you wish to use
rdkit, you should install it separately.

  .. code-block:: bash

     pip install rdkit

Updates Version 0.1.8
======================
1. The main update here is to enable mofstructure to run on Python versions
3.9 to 3.13. We have done the neccessary tests but let us know if you have
any conflicts or bugs and we will fix it.

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
