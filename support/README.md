### This sub package is for protein structure preparation only 
This is a independent package of DXGB, you don't need to conduct this under DXGB environment. 
#### Install amber tools 
```
conda install ambertools -c ambermd
```
#### Install chimera
```
conda install -c insilichem pychimera
```
#### Download and install propka 3.1
https://github.com/jensengroup/propka-3.1

#### Download and install pdb2pqr
https://sourceforge.net/projects/pdb2pqr/

#### Run test
```
pychimera prepare_structure.py
```
### Structures
input protein structure: fn + "_protein.pdb"<br>
input ligand structure: fn + "_ligand.mol2"<br>
output protein structure: fn + "_protein_prep_final.pdb"<br>
