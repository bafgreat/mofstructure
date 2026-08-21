# Changelog

All notable changes to `mofstructure` are recorded here. Versions follow the
releases published on [PyPI](https://pypi.org/project/mofstructure/).

## 0.1.9.1

### Command line

- `mofstructure_topology` accepts a directory, the way `mofstructure_database`
  and `mofstructure_curate` do. Naming a folder previously reached ASE, which
  reported it as a corrupt trajectory. A folder holding nothing readable is now
  reported as that.
- `mofstructure_topology` writes a `topology_data.csv` summary beside its
  records, pairing the two the way `porosity_data` and `fingerprint_data` are
  paired. The table is built from the merged database, not from one run.
- `status` and `source` are no longer saved to `topology_data.json`. Both
  describe the run rather than the net; both are still printed and `status`
  still sets the exit code. `mofstructure_database` drops `status` at the same
  point, so a row has one shape whichever command wrote it.
  `MOFstructure.get_topology` still returns the field.
- `mofstructure_generate_cgd` gained `--embedding`, choosing between the
  crystal's own coordinates (`deconstruction`, the default), the canonical
  barycentric placement (`ideal`), and that placement relaxed towards edges of
  equal length (`refined`). The canonical embedding carries the key, so an
  unnamed net stays identifiable from its own file.
- `mofstructure_generate_cgd` accepts `--method zeol`, which `build_cgd` has
  always implemented but the parser did not offer, and defaults to `auto`,
  reading the material from the structure. A zeolite is now named by its IZA
  framework-type code rather than by the RCSR symbol of the same net.
- Help text rewritten for both commands: worked examples, a gloss per node
  definition and per embedding, and documented exit codes. Both take
  `-v/--verbose` and `--quiet`.

### Database workflow

- `mofstructure_database` computes topology by default; `--no-topology` opts
  out. The reason it was optional, a JVM call per structure, went away when
  Systre did.
- Topology is filled in for structures already in a database rather than being
  skipped with them, matching how the fingerprint already behaved. Passing
  `-t` to an existing database previously produced nothing at all.
- `--method` accepts `auto` and defaults to it, so a folder of MOFs, COFs and
  zeolites is handled in one run instead of being forced through a MOF
  deconstruction.

### Deprecations

- `MOFstructure.get_topology` warns when `decimals` or `include_edge_centers`
  is passed. Both have been ignored since the key became exact and integral;
  they will be removed in a future release.



### Topology identification moved into Python

- Replaced the Systre call with `mofstructure.graph_net`, a pure-Python
  implementation of the Delgado-Friedrichs canonical form. **Java is no longer
  required or shipped**: the JVM jar, the `jdk4py` dependency and the
  `mofstructure.systre` module have been removed, and the wheel is about
  1.7 MB smaller.
- Validated before the switch rather than after. On the same quotient graphs,
  the Python implementation agreed with Systre on 180 of 180 curated MOFs,
  30 of 30 COFs and 10 of 10 zeolites, and with CrystalNets - an independent
  implementation with its own canonical form - on all 111 structures where
  both named a net. TD10 agrees to the digit Systre prints.
- Added `mofstructure.topology`, one call from a structure to its net. It
  accepts a file, an ASE atoms object or a CGD periodic graph, works out
  whether a framework is a MOF, a COF or a zeolite, and returns the same
  record for all three so results are comparable across chemistry.

### New

- `zeol` deconstruction for zeolites, whose vertices are the tetrahedral atoms
  and whose edges are T-O-T bridges with the oxygen contracted away. This is
  the net the zeolite literature names: ABW gives `sra`, EDI gives `edi`.
  Previously no method applied to a framework with no metal cluster to cut at.
- `hydrazide` and `ester` COF linkages, which between them cover frameworks
  that previously reported no linkage at all.
- `refine_cgd=True` on `get_topology`, and `graph_net.embedding.refined_embedding`,
  for an embedding with near-uniform edge lengths. The barycentric placement
  minimises *squared* edge length and leaves a spread of two to three between
  longest and shortest edge, which is fine to store but poor to build on;
  refinement brings `pcu` and `fcu` to 1.00. The exact embedding remains the
  default and is the reproducible one.
- `analyse_methods` reports every deconstruction that suits a material, since
  a MOF has more than one defensible net.

### Archive

- The lookup table now holds 17,454 nets, up from 2,930: the RCSR archive
  re-keyed, plus the IZA zeolite framework codes and the EPINET enumeration
  that CrystalNets distributes. Names are reported with the archive they came
  from, so an RCSR symbol is never confused with a position in a systematic
  enumeration.
- Shipped as msgpack rather than JSON, halving its size. The JSON stages are
  build intermediates and are no longer committed.

### Porosity

- **zeo++ now runs under a timeout.** A cell large enough to keep zeo++ busy
  for hours previously stalled a batch indefinitely, since the child process
  was waited on without limit. `zeo_calculation`, `get_porosity`,
  `mofstructure_porosity --timeout` and `mofstructure_database
  --porosity_timeout` all take one, defaulting to 1800 s. The child is killed
  and the structure recorded as a timeout, and its scratch files now live
  inside the directory the parent cleans up, so a killed child leaves nothing
  behind.
- **Every call returns the same keys.** A structure zeo++ cannot handle used
  to give `{}`, which the batch scripts turned into a `None` row; building the
  summary from that raised `AttributeError: 'NoneType' object has no
  attribute 'items'` *after* the whole batch had run, losing the csv. A
  failure now returns the full record with `None` in place of each number and
  a `porosity_status` of `timeout` or `failed:<code>`, so a directory of
  structures yields rows of one shape.
- **Field names are dataset ready.** `AV_Volume_fraction`, `AV_A^3`,
  `ASA_A^2`, `ASA_m^2/cm^3`, `Number_of_channels`, `LCD_A`, `PLD_A` and
  `lfpd_A` became `av_volume_fraction`, `av_a3`, `asa_a2`, `asa_m2_per_cm3`,
  `number_of_channels`, `lcd_a`, `pld_a` and `lfpd_a`. The `^` and `/` made a
  column awkward to address in pandas, SQL and parquet. **A table written by
  an earlier version uses the old names**; `porosity.LEGACY_FIELD_NAMES` maps
  one to the other. Values are plain python floats and ints rather than numpy
  scalars, which `json.dumps` refused for `Number_of_channels` without a
  custom encoder.
- zeo++ accuracy is selectable from the command line with `-a/--accuracy`,
  on `mofstructure_porosity` and `mofstructure_database`. It defaults to
  `high`, which is what both commands already used; `low` picks the cheaper
  Voronoi decomposition, worth trying before raising the timeout on a
  structure that will not finish. `mofstructure_database` also gained `-pr`,
  having previously offered no control over the porosity at all.
- `mofstructure_porosity` now honours `--probe_radius`, `--number_of_steps`
  and `--rad_file`. All three were parsed, passed to `compile_data` and then
  dropped: the call underneath took no arguments, so every run used the
  defaults whatever was asked for.
- `mofstructure_porosity` reports a structure it could not read instead of
  skipping it silently, and records it in the table with the reason.

### Cheminformatics toolkit

- **RDKit now stands in when OpenBabel is missing**, rather than the analysis
  failing. OpenBabel remains the default and nothing changes when it is
  installed. `MOFSTRUCTURE_CHEMINFO` selects a toolkit explicitly, and
  `mofdeconstructor.cheminformatics_backend()` reports which one is in use.
  Previously a missing OpenBabel printed an install hint at import and then
  raised `NameError` at the first fragment.
- Ligand names survive the switch. The IUPAC database is keyed on InChIKey,
  which the IUPAC algorithm defines rather than the toolkit, so HKUST-1 and
  MOF-5 resolve to the same names under either. Canonical SMILES is toolkit
  specific and the SMILES key simply does not match under RDKit, which is why
  it is the second key and not the first.
- The RDKit path perceives the anionic fragments deconstruction produces. A
  linker cut from its metal has no known charge, and RDKit accepts more than
  one: terephthalate is the dianion at -2 and a tetraanion with double bonds
  to the ring at -4. The least charged reading that succeeds is taken, which
  picks the chemistry. The charge RDKit reports in its own error message is
  -4 here, so that report is not used.
- A fragment RDKit cannot perceive returns empty identifiers instead of
  raising, so one awkward building unit no longer discards a deconstruction.
  The two toolkits do not always perceive a fragment identically - OpenBabel
  favours radicals where RDKit writes anions - so identifiers can differ and a
  dataset should be built with one toolkit rather than a mixture. Zeolites are
  unaffected either way, since their net comes from the T atoms and no
  cheminformatics toolkit is involved.

### Code health

- Type annotations use the built-in generics and `|` unions of PEP 585 and
  604 throughout, so `typing.Dict`, `List`, `Tuple`, `Optional` and `Union`
  are gone. The deprecated aliases warned under newer type checkers; the
  package requires python 3.10, where both forms are runtime supported.
  `Sequence`, `Iterable` and `Iterator` now come from `collections.abc`.
- Every module, and every public function and class, carries a docstring.
- The spglib call in `graph_net.symmetry` accepts both the error handling
  that returns None and the one that raises, rather than depending on which
  spglib release is installed.
- Removed dead imports, unused locals and superseded commented-out code.
  Pyflakes reports nothing across the package, tools and tests.

### Changed

- `topology_hash` now digests the canonical key rather than Systre's relaxed
  geometry. The new digest identifies the *net*: the same framework in a
  supercell, with atoms reordered or the origin moved, hashes identically,
  which the geometric digest did not. **Values stored under the old scheme
  will not match**, and `key_version` travels with the hash so a stale one is
  recognisable rather than merely wrong.
- `get_topology` returns `cgd` written from this package's own embedding. It
  describes the same net as before in a primitive rather than conventional
  cell.
- `mofstructure_topology` is now backed by the Python implementation.
  `mofstructure_systre_cgd` has been removed.
- `mofstructure_topology` writes where the rest of the package writes.
  Records are appended to `<save_dir>/Structure_Data/topology_data.json`,
  the file `mofstructure_database` already fills, so a topology run and a
  database run build one folder instead of two. `-s/--save_dir` names the
  directory, `--json` writes the full records elsewhere and `--no-save`
  prints without writing.
- **The default output directory is now `MOFstructureDB`, renamed from
  `MOFDb`.** The package identifies COFs and zeolites as readily as MOFs, so
  the old name described the output of only one third of it. Every command
  that takes `-s/--save_dir` picks the new name up from one constant,
  `filetyper.DEFAULT_SAVE_DIR`. Runs that pass `-s` explicitly are
  unaffected; runs that relied on the default will write to a new folder,
  and an existing `MOFDb` is not read or migrated - pass `-s MOFDb` to keep
  adding to it.
- `filetyper.append_json` truncates the file after rewriting it. Replacing a
  key with a shorter record previously left the tail of the old file in
  place, which made the JSON unparseable; every command that appends to a
  database file was exposed to this.

## 0.1.9.0

### Topology

- Added a chemically explicit `ligand_cluster` representation in which complete
  organic ligands and metal clusters are the two vertex classes. Edges record
  distinct periodic ligand--cluster coordination incidences; multiple donor
  bonds within one contact are consolidated without merging contacts to
  different periodic images. Polytopic and ditopic ligands remain explicit
  vertices. `collapse_ditopic=True` is available when a conventional contracted
  RCSR net is required.
- Made contracted nets invariant to atom ordering, cell origin and supercell
  choice. Contact translations for periodic rods and sheets are now reduced
  modulo the component translation lattice. This fixes order-dependent results
  for MIL-53, whose `sbus` topology now remains `pcu` across equivalent input
  representations.
- Excluded singly coordinated dangling ligands, coordinated solvent and capping
  modulators from the Systre net, preventing degree-one collisions from hiding
  the underlying framework topology.
- Made `topology_hash` independent of Systre's input-dependent node labels.
  Existing hashes from earlier releases will not match the new values. Relaxed
  coordinates may still differ by an ideal-space-group origin choice.

### Descriptors

- Added `MOFstructure.get_ligand_cluster_fingerprint()`. The fingerprint records
  ligand and cluster species, periodic connectivity, denticity, terminal
  ligands and refinement information without requiring Systre to identify the
  net. Its normalized counts and `fingerprint_hash` are invariant to atom
  ordering, cell origin and supercell expansion.

### Visualization

- Added `MOFstructure.draw_topology()` for interactive 3D inspection of the
  extracted net over the framework geometry. It supports `sbus`, `all_node`,
  `single_node` and `ligand_cluster`, periodic supercells, optional framework
  and unit-cell layers, HTML export and static image export.
- Corrected drawing across periodic boundaries by retaining graph translations,
  self-edges and distinct periodic incidences. Every displayed connection ends
  at a visible node image, and zero-length self-edges are removed.
- Added a method-specific centre-to-centre view of SBU/metal and organic/linker
  centres. Connections use solid green lines, while the abstract topology layer
  can be enabled separately with `show_topology=True`.

### Packaging and documentation

- Added the optional `draw` extra for Plotly-based visualization and documented
  the topology drawing and ligand--cluster fingerprint APIs.

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

- **Rod SBUs lost their chain connectivity.** An infinite metal-oxo rod
  (MIL-53 and similar) is periodic within itself, but that periodicity was
  discarded when the rod was contracted to a node, so the net collapsed to a
  two-dimensional `sql`.

  `method="all_node"` now splits a rod SBU into its atoms: each metal and each
  bridging carboxyl carbon becomes a node, and the oxygen atoms between them
  contract to edges, recovering the true net. MIL-53 (`Cr.cif`) gives `rna`,
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
and a `topology_data.csv` summary. At the time it was off by default because
identification then shelled out to Systre.

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
- The documented return of `get_topology()` matches what it returns. The output
  reference gained `detail`, `td10` is described as the invariant actually
  computed, and the note that `td10` and `cgd` are `None` when a keyed net has
  no ideal embedding. Where an unnamed net is concerned the README and the
  examples claimed `UNKNOWN`, which no surface returns: `get_topology()` gives
  `None` and `analyse` gives `unknown`, and both are now documented as such.
  `topology_hash` is described as the compatibility alias for `key_hash` that
  it is, below the table rather than in it, and the stored `topology_data.json`
  record is documented by the fields it actually holds.
- The Sphinx examples called `structure.MOFStructure`. The class is
  `MOFstructure`, so the snippets raised `AttributeError` when copied.

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
