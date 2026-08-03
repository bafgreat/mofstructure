<!-- markdownlint-disable MD033 MD036 -->
<div align="center">

# mofstructure

**A Python toolkit for topology, porosity and building-unit analysis of metal-organic frameworks**

[![PyPI](https://img.shields.io/pypi/v/mofstructure.svg)](https://pypi.org/project/mofstructure/)
[![Python](https://img.shields.io/pypi/pyversions/mofstructure.svg)](https://pypi.org/project/mofstructure/)
[![License: MIT](https://img.shields.io/badge/License-MIT-fca311.svg)](LICENSE)
[![Documentation](https://img.shields.io/badge/docs-github.io-14213d.svg)](https://bafgreat.github.io/mofstructure/)

<img src="https://raw.githubusercontent.com/bafgreat/mofstructure/main/docs/source/images/DUT-8.png" alt="DUT-8 framework" width="55%">

</div>
<!-- markdownlint-enable MD033 MD036 -->

`mofstructure` takes a crystal structure and answers the questions that usually
follow: what is it built from, how porous is it, and what net does it form.
It works on metal-organic frameworks, and also on covalent organic frameworks
and zeolites, from CIF or any other format ASE can read.

```python
from mofstructure import structure

mof = structure.MOFstructure(filename='UiO-66.cif')

mof.get_porosity()                 # PLD, LCD, surface area, void fraction
mof.get_sbu()                      # metal and organic secondary building units
mof.get_topology()                 # RCSR net symbol, dimensionality, TD10
mof.get_oms()                      # open metal sites
```

---

## Contents

- [Features](#features)
- [Installation](#installation)
- [Requirements](#requirements)
- [Command line](#command-line)
- [Python API](#python-api)
- [Output reference](#output-reference)
- [Documentation](#documentation)
- [Contributing](#contributing)
- [Citation](#citation)
- [License](#license)
- [Author](#author)

---

## Features

| Capability | What you get |
| --- | --- |
| **Topology** | RCSR net symbol via Systre, net dimensionality, TD10 density, a reproducible topology hash |
| **Porosity** | Pore limiting diameter, largest cavity diameter, accessible surface area and volume, channel count |
| **Guest removal** | Unbound solvent stripped automatically before every analysis |
| **Deconstruction** | Metal clusters, organic ligands, metal SBUs and organic SBUs as separate structures |
| **Cheminformatics** | SMILES, InChI and InChIKey for every building unit, plus IUPAC ligand names where known |
| **SBU characterisation** | SBU type (paddlewheel, rod-like, UiO-66-like and others) and metal coordination number |
| **Open metal sites** | Detection and local coordination environment of undercoordinated metals |
| **Periodic wrapping** | Fragments split across cell boundaries reassembled into whole molecules |

---

## Installation

From PyPI:

```bash
pip install mofstructure
```

From source, for the development version:

```bash
git clone https://github.com/bafgreat/mofstructure.git
cd mofstructure
pip install .
```

## Requirements

Python 3.10 or newer. Dependencies install automatically, with two things worth
knowing about:

- **Java** is required for topology only. Systre runs on the JVM, and the jar
  ships with the package. `mofstructure` looks for a JRE that you configure
  explicitly, then for `jdk4py`, then for `java` on your `PATH`. Every other
  feature works without it.
- **RDKit** is optional. OpenBabel handles the cheminformatics by default;
  install `rdkit` only if you want the alternative code path.

---

## Command line

Each command accepts a single structure file or a folder of them.

### Deconstruct one structure

```bash
mofstructure structure.cif                       # writes to ./MOF_building_units
mofstructure structure.cif path/to/results
```

### Build a database from a folder

```bash
mofstructure_database cif_folder                 # writes to ./MOFDb
mofstructure_database cif_folder -s path/to/results
mofstructure_database cif_folder -t              # include topology
```

Results land in `MOFDb/Structure_Data` as JSON, one file per analysis, plus a
CSV summary of the porosity. Structures already present are skipped, so an
interrupted run resumes where it stopped. Delete the output folder to force a
recomputation.

### Individual analyses

Use these when you only need one thing and want it to run fast.

```bash
mofstructure_building_units cif_folder           # deconstruction only
mofstructure_porosity cif_folder                 # porosity only
mofstructure_oms cif_folder                      # open metal sites only
```

`mofstructure_porosity` accepts a custom probe radius, cycle count and radii
file:

```bash
mofstructure_porosity cif_folder -pr 1.5 -ns 20000 -rf rad.rad
```

### Topology from the command line

```bash
mofstructure_topology structure.cif
mofstructure_topology net.cgd
mofstructure_topology ./folder
```

The net depends on how you define a node, and several definitions are
available. The same framework can legitimately give a different net under each —
for a rod MOF like MIL-53, `all_node` gives `rna` and `single_node` gives `bpq`
(matching CrystalNets), while `sbus` collapses the rod to `pcu`.

```bash
mofstructure_topology structure.cif --method all_node      # every branch point a node
mofstructure_topology ./folder --method single_node        # organic groups merged
mofstructure_topology ./folder --method sbus               # each SBU one node
mofstructure_topology ./folder --method ligand_cluster     # complete ligands and metal clusters
mofstructure_topology ./folder --method all                # all methods, one record each
```

`ligand_cluster` constructs a bipartite incidence net from the same
deconstruction as `get_ligands`, so its vertices are exactly the complete
ligands and metal clusters that deconstruction returns and every atom belongs to
one of them. An edge records coordination of a ligand to a particular periodic
image of a cluster. Multiple donor bonds to the same cluster image count as one
incidence, so chelation does not artificially increase the topological degree.

A ditopic ligand stays a vertex, which subdivides the edge it makes, and RCSR
lists no subdivided nets — UiO-66 comes back as `UNKNOWN` even though the net is
right. That is deliberate: the point of this method is how the ligands meet the
clusters, not the RCSR symbol, and the topology hash still identifies the net.
Pass `collapse_ditopic=True` to `ligand_cluster_graph` or `cgd_ligand_cluster`
to splice ditopic ligands into edges instead, which recovers the nameable net
(`fcu` for UiO-66, `pcu` for MIL-53, `tbo` either way for HKUST-1).

To ask what the ligands do rather than what the net is called, use the
fingerprint, which is read straight from the deconstruction and needs no Systre:

```python
from mofstructure import structure

mof = structure.MOFstructure(filename='UiO-66.cif')
print(mof.get_ligand_cluster_fingerprint())
```

It counts each ligand and cluster species per metal-cluster unit, with how many
clusters each ligand bridges and at what denticity, and it does not change when
the atoms are listed in another order, when the cell origin moves, or when the
same crystal is given as a supercell. That makes it sensitive to defects: a
missing linker lowers a cluster's connectivity, a linker hanging by one end is
listed under `terminal` with its own formula (which is what tells it apart from
a coordinated solvent), and a carboxylate that has dropped from bridging to
monodentate shows in the denticity histogram even though the net is unchanged.

Use `--method all` to compute every method at once. Each structure gets a single
record holding every net — nested under a `topologies` key in the JSON, and one
column group per method in the CSV — so the output drops straight into a
database:

```bash
mofstructure_topology ./folder --method all
```

For large datasets, write results to disk in batches:

```bash
mofstructure_topology ./folder --flush-every 100
```

---

## Python API

`MOFstructure` is the single entry point. Guests are removed internally, so a
structure containing solvent needs no preparation.

```python
from mofstructure import structure

mof = structure.MOFstructure(filename='UiO-66.cif')
# or pass an ASE atoms object directly
# mof = structure.MOFstructure(ase_atoms=atoms)

guest_free = mof.remove_guest()
```

### Porosity

```python
import pandas as pd

pores = mof.get_porosity(probe_radius=1.86, number_of_steps=5000, high_accuracy=True)
pd.DataFrame(pores, index=[0]).to_csv('pore.csv')
```

A structure that Zeo++ cannot analyse returns an empty dictionary rather than
raising, so a batch job is never interrupted by one difficult framework.

### Building units

```python
metal_sbus, organic_sbus = mof.get_sbu(wrap_system=True, cheminfo=True, add_dummy=False)
organic_ligands = mof.get_ligands(wrap_system=True, cheminfo=True, add_dummy=False)
```

With `cheminfo=True`, OpenBabel identifiers are attached to each fragment's
`.info` dictionary:

```python
for i, sbu in enumerate(metal_sbus):
    smi = sbu.info['smi']
    inchi = sbu.info['inchi']
    inchikey = sbu.info['inchikey']
    n_points = len(sbu.info['point_of_extension'])   # SBUs only
    sbu_type = sbu.info['sbu_type']                  # metal SBUs only
    sbu.write(f'metal_sbu_{i}.cif')
```

`add_dummy=True` marks the points of extension with dummy atoms, which makes the
cut positions explicit and easy to cap with hydrogen. Use it for SBUs only,
never when deconstructing into ligands and clusters.

### Ligand names

Building units carry identifiers but not names. To name a ligand, look it up
from its SMILES against the database that ships with the package:

```python
from mofstructure.filetyper import load_iupac_names
from mofstructure.mofdeconstructor import lookup_iupac_name

iupac_names = load_iupac_names()

_, ligands = mof.get_ligands()
for ligand in ligands:
    print(lookup_iupac_name(ligand.info['smi'], iupac_names))
# terephthalic acid
```

`lookup_iupac_name` saturates the open valences left by deconstruction and
matches on InChIKey and canonical SMILES, so the fragment does not have to be
the neutral parent molecule. It returns `None` for a ligand that is not in the
database. `mofstructure_database` does this for you and stores the result in
the `ligand_names` field of `ligands_data.json`.

### Topology from Python

```python
topology = mof.get_topology()
print(topology['topology'], topology['dimension'])
```

For finer control, drive Systre directly:

```python
from ase.io import read
from mofstructure.systre import identify_topology

identify_topology('net.cgd', input_is_cgd=True)      # from a CGD file
identify_topology('UiO-66.cif', method='all_node')   # from a structure file
identify_topology(read('UiO-66.cif'))                # from ASE atoms
```

### Drawing the net

`draw_topology` traces the net over the real framework and returns an
interactive plotly figure — nodes at the real building-unit positions, edges
following the connectivity. Needs the optional `plotly` extra
(`pip install mofstructure[draw]`).

```python
fig = mof.draw_topology(method='all_node', filename='net.html')
fig.show()   # or open net.html in a browser
```

The interactive, axis-free molecular view overlays the underlying net on the
framework and shows the unit-cell boundary. Every connection ends at a visible
centre, including connections to neighbouring periodic images. The default
green centre-to-centre network is generated by the selected topology method,
so its nodes and contractions visibly change between `sbus`, `all_node`,
`single_node` and `ligand_cluster`. Framework atoms, framework bonds and each
centre class can be toggled independently in the legend. The view shows this
method-specific network by default;
set `show_topology=True` to add the abstract blue topology edges and topology
node markers. Set `show_linker_sbu=False`,
`show_structure=False` or `show_unit_cell=False` to hide individual layers.

### Open metal sites

```python
print(mof.get_oms())
```

---

## Output reference

`get_topology()` returns:

| Key | Meaning |
| --- | --- |
| `topology` | RCSR net symbol, or `UNKNOWN` when Systre finds no match |
| `dimension` | Periodicity of the net (0, 1, 2 or 3) |
| `td10` | Topological density from Systre |
| `topology_hash` | Stable hash of the relaxed net, for indexing and duplicate detection |
| `cgd` | CRYSTAL2 text of the relaxed net |

`get_porosity()` returns:

| Key | Meaning |
| --- | --- |
| `PLD_A` | Pore limiting diameter, the largest sphere that can diffuse through |
| `LCD_A` | Largest cavity diameter, the largest sphere that fits anywhere inside |
| `lfpd_A` | Largest free sphere along the percolation path |
| `AV_A^3`, `AV_Volume_fraction` | Accessible volume and void fraction |
| `ASA_A^2`, `ASA_m^2/cm^3` | Accessible surface area |
| `Number_of_channels` | Number of distinct channels |

Custom atomic radii can be supplied through a `.rad` file, one element per line.
The extension must be `.rad` or the defaults are used silently:

```text
Mg 0.66
O 1.84
```

---

## Documentation

Full documentation is at
[docs](https://bafgreat.github.io/mofstructure/).
Release history is in [CHANGELOG.md](CHANGELOG.md).

## Contributing

Issues and pull requests are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md)
for development setup and what to include in a report.

Most problems are specific to one framework rather than general, so **please
attach the structure file** when reporting one. A CIF that reproduces the
problem is worth more than any description of it.

## Roadmap

- SBU deconstruction and topological analysis of covalent organic frameworks.

## Citation

If `mofstructure` contributes to your work, please cite:

```bibtex
@article{wonanke2026fairmofs,
  title={FAIR-MOFs: Structure-centred synthesis inference from three-dimensional
         structures of metal-organic frameworks},
  author={Wonanke, Dinga and Heine, Thomas and Longa, Antonio and others},
  year={2026},
  doi={10.21203/rs.3.rs-8375247/v1}
}
```

## License

Released under the MIT License. See [LICENSE](LICENSE).

## Author

`mofstructure` is developed by [Dinga Wonanke](https://www.dingawonanke.com).
