# deltaVinaXGB_develop
This is a machine-learning based protein-ligand scoring function.
### Setup
create environment
```
make Makefile create_environment
source activate DXG
make Makefile requirements
```
install source code
```
python setup.py install
```
Note: You still need to install rdkit(version >= 2018.03.2) and obabel 

### Usage 

Get all options

```
python run_DXGB.py --help
```
Run test

```
python run_DXGB.py --runfeatures --rw --water
```
