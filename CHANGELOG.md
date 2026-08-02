# Changelog

All notable changes to `mofstructure` are recorded here. Versions follow the
releases published on [PyPI](https://pypi.org/project/mofstructure/).

## 0.1.8.9

### Topology

Two net-construction bugs made `get_topology` report the wrong RCSR net for
whole classes of framework. Deconstruction was correct in both cases; only the
CGD handed to Systre was wrong.

- **Polytopic linkers were contracted into a clique.** Every linker was
  collapsed into edges directly between the metals it touched. For a ditopic
  linker that is one edge, which is right, but a tritopic linker (BTC and
  higher) became a triangle of edges, inflating each metal's connectivity.
  HKUST-1, a known `tbo` net, came out `reo` (its paddlewheels reading as
  8-connected instead of 4). A polytopic linker is now its own branch-point
  node joined to each metal it bridges, so HKUST-1 gives `tbo` with a net
  structurally identical to the reference (32 three-connected + 24
  four-connected vertices). Ditopic frameworks are unaffected: UiO-66 stays
  `fcu`, MOF-5 and the pillared-paddlewheel structures stay `pcu`.

- **Rod SBUs lost their chain connectivity.** An infinite metal–oxo rod
  (MIL-53 and similar) is periodic within itself, but that periodicity was
  discarded when the rod was contracted to a node, so the net collapsed to a
  two-dimensional `sql`.

  `method="all_node"` now splits a rod SBU into its atoms — each metal and each
  bridging carboxyl carbon becomes a node, and the oxygen atoms between them
  contract to edges — recovering the true net. MIL-53 (`Cr.cif`) gives `rna`,
  matching CrystalNets' AllNodes and mofid's AllNode. `method="sbus"` keeps the
  rod as a single node (giving `pcu`), so the two methods now genuinely differ
  for rods, as they should. Discrete SBUs are unaffected: HKUST-1 `tbo`,
  UiO-66 `fcu`, and pillared paddlewheels `pcu` under both methods. Validated
  against CrystalNets.

- **Added `method="single_node"`**, the coarsening CrystalNets calls
  SingleNodes. It takes the all-node net and merges each connected group of
  organic (carboxyl and linker) vertices into one vertex, leaving the metal
  vertices separate; a periodic organic group is left un-merged. MIL-53 gives
  `bpq`, again matching CrystalNets. Discrete frameworks are unchanged
  (HKUST-1 `tbo`, UiO-66 `fcu`). `mofstructure_topology` accepts it as
  `--method single_node`, and `--method all` now records sbus, all_node and
  single_node together.

### Fixes

- `find_unique_building_units(..., add_dummy=True)` crashed with an
  `IndexError`. `find_key_or_value` only matched the first atom of a broken
  bond, returning `None` for the second, and the `None` broke the ASE indexing
  that places dummy atoms. Matching is now symmetric, so a dummy is placed at
  every point of extension.
- Fixed rodlike detection for non-periodic and partially periodic structures.
- Added support for guest removal in non-periodic systems (e.g. organic cages) by retaining the heaviest connected fragment.

### Command line

- Added `--method all` to `mofstructure_topology`. It runs the node
  definitions (`sbus`, `all_node`, `single_node`) and records the nets
  in one entry per structure: the JSON nests them under a `topologies` key, and
  the CSV gives each its own columns (`sbus_topology`, `all_node_topology`,
  and so on), one row per structure, ready to load into a database.
- Made the console scripts consistent. Every script now accepts `-v/--verbose`
  (previously missing from `mofstructure_topology`, `mofstructure_systre_cgd`
  and `cof_stacking`), `mofstructure_database` accepts `--method` as the
  standard name for `--topology_method` (the old name still works), and
  `cof_stacking` accepts `-o/--output`.
- `-o` now means output everywhere. `mofstructure_database` used `-o` for
  `--oms`; that short form was removed, so use `--oms` (this is a breaking
  change for anyone who passed `-o` to that command).
- Fixed the `mofstructure_generate_cgd` entry point, which pointed at a module
  that no longer exists (`mofstructure.topology`) and failed on every install.
  It now runs. Reinstall the package to regenerate the console script.
- `cof_stacking` on a non-layered structure now prints a clear message instead
  of crashing with a traceback.
- Fixed `mofstructure_topology --finalise-only`, which required the input
  files it is meant to skip and so could never run. It now merges the existing
  batches without any inputs.

### Changes

- Removed the `ligand_cluster` topology method from `get_topology`,
  `build_cgd`, `mofstructure_topology` and `mofstructure_database`. It
  duplicated `sbus` on simple frameworks and produced non-standard (UNKNOWN)
  nets on the rest, and had no CrystalNets equivalent. Use `sbus`, `all_node`
  or `single_node`. The `ligands_and_metal_clusters` deconstruction it was
  built on is unchanged and still backs `MOFstructure.get_ligands()`. The
  `mofstructure_database` topology default is now `all_node`.
- Removed the `connect_mode` argument from `build_cgd`,
  `cgd_from_region_targets` and the `mofstructure_topology` CLI. Linker
  contraction no longer has a clique/chain choice: ditopic linkers become
  edges and polytopic linkers become nodes, unconditionally.
- Improved robustness of the deconstruction workflow and updated documentation.

## 0.1.8.7

This release fixes ligand naming, which never resolved before, and stops a
crash in the porosity code that could terminate a database run.

### IUPAC ligand names now resolve during deconstruction

`ligand_names` came back as `null` for every structure. Two separate problems
had to be fixed before a name could ever match.

The name database was scraped from PubChem, whose SMILES are generated by the
CACTVS toolkit. Canonical SMILES are only canonical within one toolkit, so a
CACTVS string never equals the OpenBabel string that `compute_smi` produces,
and a plain dictionary lookup always missed.

A ligand cut out of a framework also carries dangling valences at its points of
extension. Terephthalate leaves deconstruction as `[O]C(=O)c1ccc(cc1)C(=O)[O]`,
which is two hydrogens short of the terephthalic acid PubChem stores, so the
molecular formula, the SMILES and even the InChIKey connectivity block all
differ. Only ligands that coordinate through lone pairs, such as DABCO,
survived intact.

Both sides are now normalised the same way. Open valences are saturated with
implicit hydrogens to recover the neutral parent molecule and every molecule is
indexed under three keys: its full InChIKey, its canonical SMILES, and its
InChIKey connectivity block, which ignores protonation. The shipped database
was rekeyed accordingly.

`get_ligands()` returns ASE atoms objects, not names. Look the name up from
the SMILES on each fragment:

```python
from mofstructure import structure
from mofstructure.filetyper import load_iupac_names
from mofstructure.mofdeconstructor import lookup_iupac_name

iupac_names = load_iupac_names()

mof = structure.MOFstructure(filename='RUBTAK01.cif')
_, ligands = mof.get_ligands()

for ligand in ligands:
    print(lookup_iupac_name(ligand.info['smi'], iupac_names))
# terephthalic acid
```

The command line tools do this for you and write the result to the
`ligand_names` field of `ligands_data.json`, which previously held only
`null`.

New helpers in `mofdeconstructor`:

- `saturate_open_valences()` fills the valences left by deconstruction
- `name_lookup_keys()` builds the identifiers a molecule is indexed under
- `lookup_iupac_name()` queries the database through those identifiers

Both `mofstructure_database` and `mofstructure_building_units` use them, so
`ligand_names` is populated in the JSON output. `tools/regenerate_iupac_db.py`
rebuilds the database if it is ever refreshed from an external source.

### Canonical SMILES for building units

`compute_smi` now writes OpenBabel canonical SMILES (`can`) instead of `smi`,
which followed the input atom order. The same fragment extracted from two
different files previously produced two different strings, which made SMILES
unusable as a dictionary key.

### Porosity no longer terminates a run

zeo++ calls `abort()` when its Voronoi decomposition fails an internal volume
check. That raises SIGABRT, which no `try`/`except` can intercept, so a single
awkward framework killed the whole process. In a 114 structure test folder, 13
structures did this.

`zeo_calculation` now runs in a child interpreter and returns an empty
dictionary when the child dies, so those structures are recorded as having no
porosity data instead of ending the job. `compute_zeo_parameters` is the same
calculation without the isolation.

### Topology in the database workflow

`mofstructure_database` accepts `-t/--topology`, writing `topology_data.json`
and a `topology_data.csv` summary. It is off by default because Systre runs on
the JVM and costs a few seconds per structure.

```bash
mofstructure_database cif_folder -t
mofstructure_database cif_folder -t --topology_method all_node
```

### Other fixes

- The CSV summaries no longer raise `AttributeError` when a structure failed
  and was recorded as `null`. This affected `porosity_data.csv` as well.
- `mofstructure_generate_cgd` pointed at `mofstructure.topology:main`, a module
  that does not exist, so the command failed on every install. It now resolves
  to `mofstructure.generate_cgd:main`.
- Module level documentation across the package, and a docstring for the
  `MOFstructure` class.
- The Sphinx documentation builds without warnings, and the documented
  `get_topology()` key is now `cgd`, matching what the method returns.

## 0.1.8.6

This release introduces a major upgrade to the topology analysis workflow in `mofstructure`, providing a more robust, reproducible, and information-rich framework for topological characterization.

### Key improvements

#### 1. Enhanced topology extraction

Topology determination is now handled through a high-level interface built on top of Systre, enabling:

- Direct support for:
  - `.cgd` files
  - CIF and all ASE-readable structure formats
  - Batch processing of folders
- Automatic generation of CGD representations when needed
- Improved robustness for complex and multi-component frameworks

---

#### 2. Rich topology output

The `get_topology()` method now returns a structured dictionary containing:

- `topology` → Identified RCSR net (or `UNKNOWN`)
- `dimension` → Periodicity of the net (0D, 1D, 2D, 3D)
- `td10` → Topological density descriptor from Systre
- `topology_hash` → Stable hash of the relaxed topology
- `cgd_crystal2text` → CRYSTAL2 representation of the relaxed net

This enables reproducible identification and easy downstream storage/indexing.

---

#### 3. Relaxed-topology hashing

A deterministic topology hash is now available:

- Based on normalized relaxed coordinates
- Independent of atom ordering and numerical noise
- Suitable for:
  - database indexing
  - duplicate detection
  - large-scale screening workflows

---

#### 4. CRYSTAL2 export from relaxed topology

The topology pipeline now supports:

- Direct generation of CRYSTAL2-style CGD text from relaxed Systre output
- Optional inclusion of edge-center metadata
- Fallback conversion from original CGD when relaxed output is unavailable

---

#### 5. Memory-efficient workflow

The topology computation has been redesigned to be lightweight:

- Uses a single Systre call per structure
- Avoids redundant parsing and data duplication
- Only extracts the most informative component by default

This makes it suitable for large MOF datasets and high-throughput workflows.

---

#### 6. Improved CLI support

Topology tools now:

- Work seamlessly on files and folders
- Support CSV/JSON export of results
- Provide optional verbose output for debugging
- Maintain backward compatibility with legacy flags

---

### Example

```python
from mofstructure import structure

mof = structure.MOFstructure(filename="UiO-66.cif")
topo = mof.get_topology()

print(topo)
```

## 0.1.7

1. Implemented a robust CI/CD using git actions
2. Included add_dummy key to add dummy atoms to point of extension. This is important to effectively control the breaking point. This dummy atoms can then
   be replaced with hydrogen to fully neutralize the system.

### N.B

Be please don't use add dummy when deconstructing to ligands and clusters. The add dummy argument should be used only for sbus.
e.g

```Python
connected_components, atoms_indices_at_breaking_point, porpyrin_checker, all_regions, breaking_pairs = MOF_deconstructor.secondary_building_units(ase_atom)
metal_sbus, organic_sbus, building_unit_regions = MOF_deconstructor.find_unique_building_units(
    connected_components,
    atoms_indices_at_breaking_point,
    ase_atom,
    porpyrin_checker,
    all_regions,
    cheminfo=True,
    add_dummy=True
    )

metal_sbus[0].write('test1.xyz)
```

## 0.1.6

Added new command line tools to expedite calculations especially when working on a quite large database.

### compute only deconstruction

If you wish to only compute the deconstruction of MOFs without having to compute
their porosity and open metal sites. Then simply run the following command

```Bash
mofstructure_building_units  cif_folder
```

### compute only porosity

If you wish to only compute the porosity using default values. i.e
probe radius = 1.86, number of gcmc cycles = 10000 and default csd atomic radii, then run the following command:

```Bash
mofstructure_porosity cif_folder
```

However, if you wish to use another probe radius of maybe 1.5 and gcmc cycles of 20000 alongside custom atomic radii in a file called rad.rad, run the following command:

```Bash
mofstructure_porosity cif_folder -pr 1.5 -ns 20000 -rf rad.rad
```

### compute only open metal sites

If you are only interested in computing the open metal sites, then running the following command

```Bash
mofstructure_oms cif_folder
```

## 0.1.5

The new update enables users to include a Rad file when computing porosity using pyzeo. This allows users to specify the type of radii to use. If omitted, the default pyzeo radii will be used, which are covalent radii obtained from the CSD.

Currently, this functionality can only be used when using mofstructure as a library. This can be done as follows:

```Python
from mofstructure.porosity import zeo_calculation
from ase.io import read

ase_atom = read(filename)

pore_data = zeo_calculation(ase_atom, rad_file='rad_file_name.rad')
```

### NB

Note that filename is any ASE-readable crystal structure file, ideally a CIF file. Moreover, rad_file_name.rad is a file containing the radii of each element present in the structure file. This should be formatted as follows:

```bash
element radii
```

For example, for an MgO system, your Rad file should look like this:

```bash
Mg 0.66
O 1.84
```

Also note that of the radii file does not have the .rad extension like `rad_file_name.rad` the default radii will be used.

## 0.1.4

The new update enables the computation of open metal sites in cifs
To use this functionality run the following on the command line

```bash
mofstructure_database ciffolder --oms
```

Here ciffolder corresponse to the directory/folder containing the cif files.

After the computation the metal information will be found in a json file called `metal_info.json`. This file is found in the output folder that defaults to `MOFDb` incase none is provided.

NB

Note that computing open metal sites is computationally expensive, especially if you intend to
run it on a folder with many cif files. There I recommend that if you are not interested in computing the open metal sites simply run command without the --oms option.

```Bash
mofstructure_database ciffolder
```

This command will generate a MOFDb folder without the `metal_info.json` file. But the code will run very fast.

Also note that the `--oms` option is provided on for the `mofstructure_database` command. This is not available for `mofstructure` command which targets a single cif file. If you have a single cif file wish to compute open metal sites, simply put the cif file in a folder and rin `mofstructure_database` command on the folder (`mofstructure_database ciffolder --oms`).
