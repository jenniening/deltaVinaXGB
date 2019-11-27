### This sub package is for protein structure preparation only 
#### install amber tools 
```
conda install ambertools -c ambermd
```
#### install chimera
```
conda install -c insilichem pychimera
```
#### download and install propka 3.1
https://github.com/jensengroup/propka-3.1

#### download and install pdb2pqr
https://sourceforge.net/projects/pdb2pqr/

#### Run script
```
pychimera prepare_structure.py
```
### Structures
input protein structure: fn + "_protein.pdb"
input ligand structure: fn + "_ligand.mol2"
output protein structure: fn + "_protein_prep_final.pdb"
