import sys 
#print(sys.path)
try:
    sys.path.remove('/share/apps/chimera/1.11.2/lib/python2.7/site-packages')
except ValueError:
    pass
import os
import rdkit
import get_inputtype
from get_inputtype import get_inputtype
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import AllChem
from rdkit.ML.Cluster import Butina
import sys
if sys.platform == "linux":
    from software_path_linux import path_obabel
elif sys.platform == "darwin":
    from software_path_mac import path_obabel
import random
random.seed(10)

obable = path_obabel()


#-----------------------------------------------------------------------------
# Hyper-parameters:
# random_seed = 1
# cluster RMSD cutoff = 0.5
# global minimum energy cutoff = 1
# number conformations = 1000
#-----------------------------------------------------------------------------



def gen_conformers(mol, numConfs):
    """
    Generate conformation using ETKDG
    """

    
    ps = AllChem.ETKDG()
    ps.randomSeed = 1
    ps.pruneRmsThresh = 0.1
    ps.numThreads = 0
    ids = AllChem.EmbedMultipleConfs(mol,numConfs, ps)
    print(list(ids))
    
    # minimization
    for conformerId in ids:
        mp = AllChem.MMFFGetMoleculeProperties(mol)
        # solvent(water) dielectronic constant
        mp.SetMMFFDielectricConstant(80)
    
        ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=conformerId)
        ff.Initialize()
        ff.Minimize(maxIts=1000)
    
    return list(ids)

def calc_energy(mol, conformerId):
    """
    MMFF minimization and energy calculation 
    """
    mp = AllChem.MMFFGetMoleculeProperties(mol)
    # solvent(water) dielectronic constant
    mp.SetMMFFDielectricConstant(80)
    
    ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=conformerId)
    energy = ff.CalcEnergy()
    
    return energy

def cluster_conformers(mol, mode="RMSD", threshold=0.2):
    """
    Cluster conf based on heavy atom rmsd and Butina is used for clustering
    """

    # get heavy atom idx
    heavyatomidx = []
    for a in mol.GetAtoms():
            if a.GetAtomicNum() != 1:
                heavyatomidx.append(a.GetIdx())
    heavyatomidx = set(heavyatomidx)

    # align on heavy atom for each pair and get dmat
    n = mol.GetNumConformers()
    print("The total number of conformers before clustering: " + str(n))

    dmat = []
    for i in range(n-1):
        for j in range(i+1,n):
            dmat.append(Chem.rdMolAlign.AlignMol(mol,mol,i,j,atomMap = [(k, k) for k in heavyatomidx]))
            #dmat.append(Chem.rdMolAlign.GetBestRMS(mol,mol,i,j,map=[[(k, k) for k in heavyatomidx]]))
            #print(dmat)
    # clustering
    rms_clusters = Butina.ClusterData(dmat,n,
                                      threshold, isDistData=True, reordering=True)
    return rms_clusters


def runGenerator(fn, input_file, outfile, outfile_min, numConfs,addH = False):
    """
    Generate conformation as sdf for each mol2 file
    """

    ftype = get_inputtype(input_file)
    
    w = Chem.SDWriter(outfile)
    w_min = Chem.SDWriter(outfile_min)
    
    if ftype == "mol2":
        mol = Chem.MolFromMol2File(input_file,removeHs=False)
    elif ftype == "smi":
        mol = Chem.MolFromSmiles(input_file,removeHs=False)
    elif ftype == "sdf":
        mol = Chem.SDMolSupplier(input_file,removeHs=False)[0]

    if addH:
        mol = Chem.AddHs(mol)
    print("mol has been read")


    # generate the confomers
    conformerIds = gen_conformers(mol, numConfs)
    print(conformerIds)
    # clustering
    rmsClusters = cluster_conformers(mol, "RMSD", 0.5)
    
    conformerPropsDict = {}
    
    n = 0
    energy_list = []
    for clusterId in rmsClusters:
        n = n + 1
        for conformerId in clusterId[:1]:
            conformerPropsDict[conformerId] = {}
            # energy minimise (optional) and energy calculation
            conformerPropsDict[conformerId]["energy_abs"] = calc_energy(mol, conformerId)
            energy_list.append(conformerPropsDict[conformerId]["energy_abs"])
            conformerPropsDict[conformerId]["SMILE"] = Chem.MolToSmiles(mol)
            conformerPropsDict[conformerId]["cluster_no"] = n
            conformerPropsDict[conformerId]["pdb_id"] = fn
            conformerPropsDict[conformerId]["initial_id"] = conformerId
            for key in conformerPropsDict[conformerId].keys():
                mol.SetProp(key, str(conformerPropsDict[conformerId][key]))
            w.write(mol, confId = conformerId)
    
    print("The total number of conformers after clustring: " + str(n))
    print("The lowest energy: " + str(min(energy_list)))
    
    # get conformation with lowest energy
    for clusterId in rmsClusters:
        n += 1
        for conformerId in clusterId[:1]:
            if conformerPropsDict[conformerId]["energy_abs"] == min(energy_list):
                for key in conformerPropsDict[conformerId].keys():
                    mol.SetProp(key, str(conformerPropsDict[conformerId][key]))

                w_min.write(mol, confId = conformerId)
    w.flush()
    w.close()
    w_min.flush()
    w_min.close()

    return None

def minimize_native(native):
    """
    Get the local minimum structure
    :return: local minimum structure
    """
    file_type = get_inputtype(native)
    m = None
    if file_type == "mol2":
        m = Chem.MolFromMol2File(native,removeHs=False)
    elif file_type == "sdf":
        m = Chem.SDMolSupplier(native,removeHs=False)[0]
    else:
        print("Error: Native File Type Is Wrong !")
    if m == None:
        print("Error: Native File Type Is Wrong !")
        
    mp = AllChem.MMFFGetMoleculeProperties(m)
    # solvent(water) dielectronic constant
    mp.SetMMFFDielectricConstant(80)
    
    ff = AllChem.MMFFGetMoleculeForceField(m, mp)
    ff.Initialize()
    ff.Minimize(maxIts=1000)
    
    return m,mp

def get_lowest_energy(lowest):
    """
    Get the global minimum energy
    """
    df_confs = PandasTools.LoadSDF(lowest)
    df_confs["energy_abs"] = df_confs["energy_abs"].astype(float)
    lowest = df_confs.sort_values(["energy_abs"]).energy_abs.min()

    return lowest


def get_native_energy(native, write = True):
    """
    Get the local minimum energy
    """
    m,mp = minimize_native(native)
    ff = AllChem.MMFFGetMoleculeForceField(m,mp)
    native_energy = ff.CalcEnergy()
    
    # write local_min of native
    if write:
        w = Chem.SDWriter(native.split(".")[0] + "_local_min.sdf")
        conformerPropsDict = {}
        # energy minimise (optional) and energy calculation
        conformerPropsDict["energy_abs"] = native_energy
        conformerPropsDict["SMILE"] = Chem.MolToSmiles(m)
        for key in conformerPropsDict.keys():
            m.SetProp(key, str(conformerPropsDict[key]))
        w.write(m)
        w.flush()
        w.close() 

    return m, native_energy

def energy_difference(lowest, native_energy):
    """
    Get the dE feature
    :param lowest: global minimum energy
    :param native_energy: local minimum energy
    """
    dE = lowest - native_energy

    return dE

def num_structure_change(confs,native):
    """
    Get number of conformations satisfying requirements --> for entropy
    """
    df_confs = PandasTools.LoadSDF(confs)
    df_confs["energy_abs"] = df_confs["energy_abs"].astype(float)
    lowest = df_confs.sort_values(["energy_abs"]).energy_abs.min()
    num_1 = df_confs[df_confs["energy_abs"] < lowest + 1.0].shape[0]
    num_2 = df_confs[df_confs["energy_abs"] < native].shape[0]
    
    return num_1, num_2

def get_RMSD(local_min,lowest,fn):
    """
    Get the RMSD between local minimum and global minimum
    """

    RMSD = None
    lowest = Chem.SDMolSupplier(lowest)[0]

    heavyatomidx = []
    for a in lowest.GetAtoms():
            if a.GetAtomicNum() != 1:
                heavyatomidx.append(a.GetIdx())
    
    m = Chem.RemoveHs(local_min)
    try:
        RMSD = Chem.rdMolAlign.GetBestRMS(m,lowest)
    except:

        RMSD = Chem.rdMolAlign.AlignMol(m,lowest,-1,-1,atomMap = [(k, k) for k in heavyatomidx])
        #out_error.write(fn + "\n")
    return RMSD

def feature_cal(outfile,fn, native, datadir, calc_type = "GenConfs"):
    confs = os.path.join(datadir, fn + "_ligand_confs.sdf")
    lowest = os.path.join(datadir, fn + "_ligand_global_min.sdf")
    if calc_type == "GenConfs":
        numConfs = 1000
        if not (os.path.exists(confs) and os.path.exists(lowest)):
            runGenerator(fn, native, confs, lowest, numConfs)

    lowest_energy = get_lowest_energy(lowest)
    local_min,native_energy = get_native_energy(native)
    dE = energy_difference(lowest_energy, native_energy)
    RMSD = get_RMSD(local_min,lowest,fn)
    num_1, num_2 = num_structure_change(confs,native_energy)
    outfile.write(fn + "," + str(dE) + "," + str(RMSD) + "," + str(num_1) + "," + str(num_2) + "\n")
    outfile.close()


    
    return None


def main():
    args = sys.argv[1:]
    if args[-1] == "file":
        pdbfile = open('%s'%(sys.argv[1] + sys.argv[2]),'r')
        pdblist = []
        for i in pdbfile.readlines():
            pdblist.append(i[0:4])
    else:
        pdblist = []
        pdblist.append(sys.argv[2])
    datadir = sys.argv[1]
    print(datadir)
    outfile = open(datadir + "dE_RMSD.csv","w")
    outfile.write("pdb,dE_global,RMSD_global,number0,number1\n")
    for fn in pdblist:
        if fn + "_ligand.sdf" in os.listdir(datadir):
            inlig = os.path.join(datadir,fn + "_ligand.sdf")
            try:
                feature_cal(outfile,fn,inlig,datadir)
            except:
                print("Ligand structure should be checked")
        elif fn + "_ligand.mol2" in os.listdir(datadir):
            inlig =  datadir + fn + "_ligand.mol2"
            try:
                feature_cal(outfile,fn,inlig,datadir)
            except:
                print("Ligand structure should be checked")
        elif fn + "_ligand_rigid.pdbqt" in os.listdir(datadir):
            # Example: convert pdbqt type into sdf, bonding infor should be checked
            olddir = os.getcwd()
            os.chdir(datadir)
            inlig =  datadir + fn + "_ligand.pdbqt"
            inlig_out =  datadir + fn + "_ligand.sdf"
            os.system(obable + " -ipdbqt " + datadir + inlig + " -osdf -O " + datadir + inlig_out)
            os.chdir(olddir)
            try:
                feature_cal(outfile,fn,inlig_out,datadir)
            except:
                print("Ligand structure should be checked")
        else:
            print("wrong ligand input file format, it should be mol2, sdf, pdbqt")
    outfile.close()

    return None



if __name__ == "__main__":

    main()
