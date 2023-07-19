# Generalities 
This is an elaborate python module that provides simple fucntions for
manipulation metal organic frameworks and  other porous systems such as 
COFs and Zeolites. The module can enable 

1. Computation of geometry properties of MOFs. It calls zeo++ in the background and enables a quick computation of all porosity information and returns a python dictionary object.

2. Automated removal of unbound guest molecules

3. Deconstruction of metal organic frameworks into building units. And for each buidling units (organic ligand, metal cluster, organic sbu and metal sbu) computes their cheminformatic identifiers such as SMILES strings, inchi and inchikey 

4. wraps systems around unit cell so as to remove effect of pbc. 

5. Seperation of building units into regions 

# Installation 
## Option 1 
 `pip install mofstructure`
## Option 2 
```
  git clone https://github.com/bafgreat/mofstructure.git mofstructure
  cd mofstructure
  pip install -e .
```
# How to: Quick start 
## Run on the commandline
### Building units 
Simply run the following command on a cif file or any ase readable file format containing a MOF. 
```
mofstructure cif_file 
```
The script will deconstruct the MOF present in the cif file and load the output in a folder called 'MOF_building_units' in the current directory/folder. If you wish to load the output in a specific folder, simply add the path to the folder as follows:
```
mofstructure cif_file path_to_result_folder
```
For multiple cif files. Simply run a loop and all the Results will be saved in the 
```
for cifs in ciffiles:
    mofstructure cifs path_to_result
```
### Creating a database
If you have a folder containg many cif files for different MOF, you could easily create a database. To create such a database, simply run the following command. 
```
mofstructure_database ciffolder
```
Here the 'ciffolder' is the folder containing the cif files. The ouput will be saved in the default folder called 'MOFDb' in the current folder. Again you can choose the path to the save folder by simply listing it at the end of the command. 
```
mofstructure_database ciffolder path_to_result
```
## Use as a libray 
1. Importing module

```
from  mofstructure import mofdeconstructor
from  mofstructure import porosity
from  mofstructure import buildingunits
from ase.io import read, write 
import pandas as pd
```
2. Read using ase 
```
ase_atom = read(cif_file)
```
3. Removal of unbound guest 
```
no_guest_indices = mofdeconstructor.remove_unbound_guest(ase_atom)
no_guest_atom = ase_atom[no_guest_indices]
```
4. Computing porosity 
```
pores = porosity.zeo_calculation(ase_ato, probe_radius=1.86, number_of_steps=5000)
df = pd.DataFrame(pores, index=[0])
df.to_csv('pore.csv')
```
5. sbus and ligands 
```
connected_components, atoms_indices_at_breaking_point, porpyrin_checker, all_regions = MOF_deconstructor.secondary_building_units(ase_atom)

metal_sbus, organic_sbus, building_unit_regions = MOF_deconstructor.find_unique_building_units(
    connected_components,
    atoms_indices_at_breaking_point,
    ase_atom, 
    porpyrin_checker,
    all_regions,
    cheminfo=True
    )
```
### when cheminfo = True
openbabel is called to compute all chemifomatic information,
which are all stored on the ase_atom.info
metal_sbus and  organic_sbus  list that contains all the unique instances of the metal sbus and organic sbus. 
### extracting cheminfor 
For each instance in a building unit the various chemiformatic informations are as follows. 
```
for i,  sbu in enumerate(metal_sbu):
    smi = sbu.info['smi]
    inchi = sbu.info['inchi']
    inchikey = sbu.info['inchikey]
    # for sbus only 
    number_of_point_of_extension = sbu.info['point_of_extension']
    #for metal sbus only 
    sbu_type = sbu.info['sbu_type'] # sbu_type :rodlike, irmof, uoi66, paddlewheel e.t.c
    # write 
    sbu.write('metal_sbu_'+str(i)+'.cif')
```
# Summary 
All of above codes can be run by a single function found mofstructure.buildingunits
## for a single cif file
``` 
from  mofstructure import buildingunits
buildingunits.work_flow(ciffile, result_folder)
```
## for multiple cif files 
```
from  mofstructure import buildingunits
import glob 
all_cifs = glob.glob(folder/*cif)
# folder corresponds to the folder containing all the cif files
for cif_files in all_cifs:
    buildingunits.work_flow(cif_files, result_folder)
# Note that result_folder can be any path. If the path does not exist, it will create one and populate it with all the data. 
```
# Support 
The module contains much more functionalities. If you are struggling or if you wish to compute a new quantity that is not yet present, feel free to send me an email. dak52@uclive.ac.nz 

# Roadmap
In the future the code should be able to: 
1. Compute rcsr topological code
2. Subsitutue building units in a MOF to enable framework functionalisation
3. Automatic curation of cifs
4. Decontsruction of COFs into their building units 


