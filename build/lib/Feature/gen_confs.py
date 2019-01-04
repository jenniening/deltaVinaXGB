"""
Generate Ligand Conformations by RDkit
"""
__author__ = "Jianing Lu"
__copyright__ = "Copyright 2018, NYU"
__license__ = ""

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

from rdkit import Chem
from rdkit.Chem import AllChem, TorsionFingerprints
from rdkit.ML.Cluster import Butina
from get_inputtype import get_inputtype

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------


def gen_conformers(mol, numConfs=100, maxAttempts=1000,
                   pruneRmsThresh=0.1,
                   useExpTorsionAnglePrefs=True,
                   useBasicKnowledge=True,
                   enforceChirality=True):
    """Generate conformation with MMFF opt"""
    # generate conf
    ids = AllChem.EmbedMultipleConfs(mol, numConfs=numConfs,
                                     maxAttempts=maxAttempts,
                                     pruneRmsThresh=pruneRmsThresh,
                                     useExpTorsionAnglePrefs=useExpTorsionAnglePrefs,
                                     useBasicKnowledge=useBasicKnowledge,
                                     enforceChirality=enforceChirality,
                                     numThreads=0)
    # MMFF optimize
    for cid in ids:
        _ = AllChem.MMFFOptimizeMolecule(mol, confId = cid)

    return list(ids)

def calc_energy(mol, conformerId):
    """Calculate MMFF energy"""

    mp = AllChem.MMFFGetMoleculeProperties(mol)
    ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=conformerId)
    results = {}
    results["energy_abs"] = ff.CalcEnergy()

    return results


def runGenerator(fn,input_file, output_file ,numConfs, maxAttempts, pruneRmsThresh, addH = False):
    """Generate conformation as sdf for each mol2 file"""
    type = get_inputtype(input_file)
    w = Chem.SDWriter(output_file)
    if type == "mol2":
        m = Chem.MolFromMol2File(input_file,removeHs=False)
    elif type == "smi":
        m = Chem.MolFromSmiles(input_file,removeHs=False)
    elif type == "sdf":
        mol = Chem.SDMolSupplier(input_file,removeHs=False)
        for m1 in mol:
            m = m1


    if addH:
        m = Chem.AddHs(m)


    # generate the confomers
    conformerIds = gen_conformers(m, numConfs, maxAttempts,
                                      pruneRmsThresh, True, True, True)
    conformerPropsDict = {}
    for conformerId in conformerIds:
        # energy minimise (optional) and energy calculation
        conformerPropsDict[conformerId]= calc_energy(m, conformerId)
        conformerPropsDict[conformerId]["SMILE"] = Chem.MolToSmiles(m)
        conformerPropsDict[conformerId]["pdb_id"] = fn
        for key in conformerPropsDict[conformerId].keys():
                m.SetProp(key, str(conformerPropsDict[conformerId][key]))
        w.write(m, confId=conformerId)
    w.flush()
    w.close()

    return None

def RMSD(initial,m):
    """ Get RMSD """
    m1 = Chem.RemoveHs(initial)
    m2 = Chem.RemoveHs(m)
    RMS = AllChem.GetBestRMS(m1, m2)

    return RMS

def cluster_conformers(mol, mode="RMSD", threshold=0.2):
    """Cluster conf based on heavy atom rmsd and Butina is used for clustering"""

    # get heavy atom idx
    heavyatomidx = []
    for m in mol:
        for a in m.GetAtoms():
            if a.GetAtomicNum() != 1:
                heavyatomidx.append(a.GetIdx())
    heavyatomidx = set(heavyatomidx)

    # align on heavy atom for each pair and get dmat
    n = 0
    for num, m in enumerate(mol):
        n = n + 1
    print("The total number of conformers before clustering: " + str(n))

    dmat = []
    for i in range(n-1):
        for j in range(i+1,n):
            dmat.append(Chem.rdMolAlign.AlignMol(mol[i],mol[j],atomMap = [(k, k) for k in heavyatomidx]))
    # clustering
    rms_clusters = Butina.ClusterData(dmat,n,
                                      threshold, isDistData=True, reordering=True)
    return rms_clusters

def energy(mol):
    """ Get Energy """
    mp = AllChem.MMFFGetMoleculeProperties(mol)
    ff = AllChem.MMFFGetMoleculeForceField(mol,mp)
    energy = ff.CalcEnergy()

    return energy

def Get_infor(fn, mol, outfile_final, outfile_min):
    """ Get RMSD_native, RMSD_min, dE_native, dE_min infor """
    w = Chem.SDWriter(outfile_final)
    w_min = Chem.SDWriter(outfile_min)
    # clustering
    rmsClusters = cluster_conformers(mol, "RMSD", 0.2)
    conformerPropsDict = {}
    n = 0
    energy_list = []
    for clusterId in rmsClusters:
        n = n + 1
        for conformerId in clusterId[:1]:
            conformerPropsDict[conformerId] = {}
            conformerPropsDict[conformerId]["energy_abs"]= energy(mol[conformerId])
            energy_list.append(energy(mol[conformerId]))
            conformerPropsDict[conformerId]["SMILE"] = Chem.MolToSmiles(mol[conformerId])
            conformerPropsDict[conformerId]["cluster_no"] = n
            conformerPropsDict[conformerId]["pdb_id"] = fn
            conformerPropsDict[conformerId]["initial_id"] = conformerId
            # only can add properties after copy the initial one
            m = mol[conformerId]
            for key in conformerPropsDict[conformerId].keys():
                m.SetProp(key, str(conformerPropsDict[conformerId][key]))
            w.write(m)
    print("The total number of conformers after clustring: " + str(n))
    print("The lowest energy: " + str(min(energy_list)))
    for clusterId in rmsClusters:
        n += 1
        for conformerId in clusterId[:1]:
            if conformerPropsDict[conformerId]["energy_abs"] == min(energy_list):
                m = mol[conformerId]
                for key in conformerPropsDict[conformerId].keys():
                    m.SetProp(key, str(conformerPropsDict[conformerId][key]))

                w_min.write(m)

    w.flush()
    w.close()
    w_min.flush()
    w_min.close()


def run_gen_conformations(fn,inlig, datadir):
    outfile = datadir + inlig.split(".")[0] + "_conformers.sdf"
    infile = datadir + inlig
    runGenerator(fn,infile,outfile,300,1000,0.1)
    conformer_file = outfile
    conformers = Chem.SDMolSupplier(conformer_file, removeHs=False)
    outfile_final = datadir + inlig.split(".")[0] + "_conformers_final.sdf"
    outfile_min = datadir +  inlig.split(".")[0]  + "_conformers_min.sdf"
    Get_infor(fn,conformers, outfile_final,outfile_min)

if __name__ == "__main__":

    testdir = '../Test/'
    inlig = '3ao4_ligand.mol2'
    fn = "3ao4"
    run_gen_conformations(fn,inlig,testdir)
