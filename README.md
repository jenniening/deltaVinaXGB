# deltaVinaXGB
This is a machine-learning based protein-ligand scoring function.
### Setup
Create environment
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
Note: if pandas/numpy can't be imported after installing rdkit, just run the make Makefile requirements command again. <br>

To calculate Vina and SASA features, you should install mgltools, msms and a modified version of AutoDock Vina. To obatin deltaVinaRF predicted scores,you should also install R and its randomForest library. Detailed information can be found in our deltaVinaRF Tutorial http://www.nyu.edu/projects/yzhang/DeltaVina/tutorial.html <br>
Remember to change the atmtypenumbers file in msms into our provided atmtypenumbers file (in current directory). <br>

Install source code
```
python setup.py install
```
Before running model, remember change the python, R, obabel, mgl tools paths into your own directorys in Feature/software_path_mac.py or Feature/software_path_linux.py file. <br>

### Prepare Data
Before calculating features, three structure inputfiles are needed:<br>
pdbid_ligand.mol2/sdf         --> ligand structure file<br>
pdbid_protein.pdb             --> protein structure file<br>
pdbid_protein_all.pdb         --> protein with water molecules structure file<br>

Note: these three files are needed to run predictions. All files should include hydrogens. If protein files with water molecules are not available, just copy the original protein.pdb to protein_all.pdb. 

If features have been already calculated, only features files are needed: <br>
Input.csv                     --> Input feature file <br>


### Run model

All scripts are in Feature directory.<br>

```
cd Feature
```
Check all options and defaults

```
python run_DXGB.py --help
```
Predict scores for Input.csv
```
python run_DXGB.py --runrf --average
```
--runrf is for deltaVinaRF scores <br>
--average is for ensemble predictions from 10 models.<br>

Calculate features and predict scores
```
python run_DXGB.py --runfeatures --datadir ../Test_2al5 --pdbid 2al5 --average
```
--runfeatures is for feature calculation, default is to calculate all features.<br>
--datadir is for structure files datadir.<br>
--pdbid is for structure pdbid, can be other types of index.<br>
Default is only calculating scores using original ligand stucture. <br>
If you want to get structures with water molecules, and optimized ligands
```
python run_DXGB.py --runfeatures --datadir ../Test_2al5 --pdbid 2al5 --water rbw --opt rbwo --average 
``` 
--water is for consideration of water, rbw is to consider both receptor-bound water and bridging water. <br>
--opt is for optimization, rbwo is to optimize ligand in no water, bridging water, and receptor-bound water environments.<br>

The predicted scores for different structures of Vina, and deltaVinaXGB will be saved in outfile (default is score.csv) in datadir.<br>
If you want to get deltaVinaRF scores as well, add --runrf. <br>

Note:
1) Ligand structure should includes both atom and bond information, such as mol2 and sdf. Be careful when using mol2 file as input format, some atom types are not recognized in RDKit (O.co2 for O in C-PO32- group). 
2) Using different version of RDKit, the ligand stability features can be slightly different.
3) Abbrevations: RW --> receptor water; BW --> bridging water

### Reference
1. Wang, C.; Zhang, Y. K., Improving Scoring-Docking-Screening Powers of Protein-Ligand Scoring Functions using Random Forest. J. Comput. Chem. 2017, 38, 169-177. https://doi.org/10.1002/jcc.24667
2. deltaVinaRF http://www.nyu.edu/projects/yzhang/DeltaVina/




