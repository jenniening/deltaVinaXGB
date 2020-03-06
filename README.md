# deltaVinaXGB
This is a machine-learning based protein-ligand scoring function.
### Setup
Create environment
```
make Makefile create_environment
conda activate DXGB
make Makefile requirements
```
Install xgboost, rdkit and obabel, these packages can be easily installed using anaconda
```
conda install -c conda-forge xgboost=0.80.0
conda install -c rdkit rdkit=2019.03.1
conda install -c openbabel openbabel
```
Install source code
```
python setup.py install
```
### All Dependencies 

To calculate Vina and SASA features, you should install mgltools, msms and a modified version of AutoDock Vina. To obatin deltaVinaRF predicted scores, you should also install R and its randomForest library. <br>

Download MGLTools and MSMS from http://mgltools.scripps.edu/downloads, choose right version for your platform (Linux or Mac).<br>
#### Install MGLTools
```
tar -xvzf mgltools_x86_64Linux2_1.5.6.tar.gz 
cd mgltools_x86_64Linux2_1.5.6/ 
./install.sh
```
Remember use correct mgltoos name (Linux or Mac) when install mgltools.
#### Install msms
```
mkdir msms 
tar -xvzf msms_i86_64Linux2_2.6.1.tar.gz -C msms 
cd msms 
cp msms.x86_64Linux2.2.6.1 msms 
```
In msms folder, there is a script pdb_to_xyzr. Change the line numfile = "./atmtypenumbers" to be numfile = "YourPATHofddeltaVinaXGB/DXGB/atmtypenumbers". atmtypenumbers file we used can be found in deltaVinaXGB/DXGB directory <br><br>
Test pdb_to_xyzr
```
pdb_to_xyzr 1crn.pdb > 1crn.xyzr
```
If it doesn't work, try 
```
./pdb_to_xyzr 1crn.pdb > 1crn.xyzr
```
If error (nawk: command not found) appears, change nawk to awk in pdb_to_xyzr (line 31) <br>

#### Install modified AutoDock Vina
Download modified Vina from  https://github.com/chengwang88/vina4dv <br>
Extract files from zip file. <br>
If directory name is vina4dv_master, change it into vina4dv. <br>

#### Install R (Only needed for deltaVinaRF)
Download R from https://cran.r-project.org/ and install randomForest in R by
```
install.packages('randomForest')
```
#### Set the environment variable
If you have the dependencies installed already. Several environment variables need to be set in .bashrc (Linux) or .bash_profile (macOS) file in your home directory. An example is given below. You can modify the path based on your case. In this example, all softwares are installed under /home/jl7003 directory.<br>
```
# path for MSMS 
export PATH=$PATH:/home/jl7003/msms/

# set mgltool variable (if mac, should change mgltools_x86_64Linux2_1.5.6 into your downloaded mac version)
export PATH=$PATH:/home/jl7003/mgltools_x86_64Linux2_1.5.6/bin/
export MGL=/home/jl7003/mgltools_x86_64Linux2_1.5.6/ 
export MGLPY=$MGL/bin/python 
export MGLUTIL=$MGL/MGLToolsPckgs/AutoDockTools/Utilities24/ 

# set vina dir  (if mac, should use /mac/release/)
export VINADIR=/home/jl7003/vina4dv/build/linux/release/ 

# set DXGB dir (needed if run deltaVinaRF20 score)
export DXGB=/home/jl7003/deltaVinaXGB/DXGB
```
After set all environment variables, open a new window to enable all setups. <br>

### Prepare Data
Before calculating features, three structure inputfiles are needed:<br>
pdbid_ligand.mol2/sdf         --> ligand structure file<br>
pdbid_ligand.pdb              --> (optional, needed to provid if you are not sure about your mol2/sdf file quality; this pdb file can be used to calculate all features excluding ligand stability features, and you still can get score) ligand structure file <br>
pdbid_protein.pdb             --> protein structure file<br>
pdbid_protein_all.pdb         --> protein with water molecules structure file<br>
Examples for structure files can be found in Test_2al5 directory.<br>
Note: these three files are needed to run predictions. All files should include hydrogens. If protein files with water molecules are not available, just copy the original protein.pdb to protein_all.pdb. 

If features have been already calculated, only features files are needed: <br>
Input.csv                     --> Input feature file <br>
Example for Input.csv can be found in Test directory.<br>

### Run model
Activate DXGB environment<br>
```
conda activate DXGB
```
All scripts are in DXGB directory.<br>
```
cd DXGB
```
Check all options and defaults<br>
```
python run_DXGB.py --help
```
The script can be run for one complex by
```
python run_DXGB.py --runfeatures --datadir ../Test_2al5 --pdbid 2al5 --average
```
--runfeatures is feature calculation, default is to calculate all features<br>
--datadir is for structure files datadir<br>
--pdbid is for structure pdbid, can be other type of index<br>
--average is to calculate average scores from 10 models<br>
Or it can also be run by providing a list of protein-ligand complex with input features as in Input.csv
```
python run_DXGB.py --datadir ../Test --average 
```
Default is to predict scores for provided structures. If you want to get scores with explicit water molecules, and optimized ligands:
```
python run_DXGB.py --runfeatures --datadir ../Test_2al5 --pdbid 2al5 --water rbw --opt rbwo --average 
``` 
--water is for consideration of water, **rbw** is to consider both receptor-bound water and bridging water. Other options: rw: only consider receptor-bound water; bw: only consider bridging water, using bw<br>
--opt is for optimization, **rbwo** is to optimize ligand in no water, bridging water, and receptor-bound water environments.
Other options: o: optimize ligand without water; rwo: optimize ligand with no water and with receptor-bound water; bwo: optimize ligand with no water and bridging water. Reminder optimization is based on protein structure, thus, you should optimization option that matches your water option.<br>

The calculated features will be saved in <code>Input.csv</code> file<br>
The predicted scores for different structures of Vina, and deltaVinaXGB will be saved in outfile (default is <code>score.csv</code>) in datadir.<br>
If you want to get deltaVinaRF scores as well, add --runrf. <br>

Note:
1) Linux is the recommended system for our package. Our package also supports Mac OSX. 
2) Ligand structure should includes both atom and bond information, such as mol2 and sdf. Be careful when using mol2 file as input format, some atom types are not recognized in RDKit (O.co2 for O in C-PO32- group). 
3) Using different version of RDKit, the ligand stability features can be slightly different.
4) If the mol2/sdf file can not be processed by RDKit, we will skip the ligand stability feature calculation and used default value (dE:-300, RMSD:300). 
4) Abbrevations: RW --> receptor water; BW --> bridging water

### Reference
1. Lu, J. N.; Hou, X. B.; Wang, C.; Zhang, Y. K., Incorporating Explicit Water Molecules and Ligand Conformation Stability in Machine-Learning Scoring Functions. J. Chem. Inf. Model. 2019, https://pubs.acs.org/doi/abs/10.1021/acs.jcim.9b00645
2. Wang, C.; Zhang, Y. K., Improving Scoring-Docking-Screening Powers of Protein-Ligand Scoring Functions using Random Forest. J. Comput. Chem. 2017, 38, 169-177. https://doi.org/10.1002/jcc.24667





