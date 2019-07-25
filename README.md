# deltaVinaXGB_develop
This is a machine-learning based protein-ligand scoring function.
### Setup
create environment
```
make Makefile create_environment
source activate DXG
make Makefile requirements
```
You still need to install rdkit(version >= 2018.03.2) and obabel, these two packages can be easily installed using anaconda

```
conda install -c rdkit rdkit
conda install -c openbabel openbabel
```
Note: if pandas/numpy can't be imported after installing rdkit, just rerun the make Makefile requirements command. 

install source code
```
python setup.py install
```

### Usage 

Run model

Check all options 
```
python run_DXGB.py --help
```
Run test example

```
python run_DXGB.py --runfeatures --average
```
--runfeatures flag is for feature calculation, default is to calculate all features.
--average flag is for average predictions from 10 models.

