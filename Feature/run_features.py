import os
import sys

import rdkit 
from rdkit import Chem
import get_pdbinfo

import CyclicPoseFragmentation_noH
from CyclicPoseFragmentation_noH import runMethod
import get_fragment_Vina
from get_fragment_Vina import generate_data
from get_fragment_Vina import run_Vina_Fragment

import bw
from bw import cal_BW
import rw
from rw import get_Crw
import opt
from opt import get_Co
import cal_vina58
from cal_vina58 import featureVina
import cal_sasa
from cal_sasa import cal_SASA
import cal_dERMSD
from cal_dERMSD import feature_cal
import cal_ion
from cal_ion import cal_Ni

import combine_data
from combine_data import combine


if sys.platform == "linux":
    from software_path_linux import path_obabel
elif sys.platform == "darwin":
    from software_path_mac import path_obabel

obabel = path_obabel()

def run_fragments(fn, datadir, inlig, inlig_pdb, water = True, decoy = False):
    '''
    get Vina58 for core and side

    '''
    ### decoy = True only for structures have many docking poses ###
    olddir = os.getcwd()
    os.chdir(datadir)
    cmd = "mkdir Frags"
    os.system(cmd)
    os.chdir("Frags")

    ### run fragmentation ###
    infmt = inlig.split(".")[1]
    if "/" in inlig:
        basename = inlig.split("/")[-1].split(".")[0]
    else:
        basename = inlig.split(".")[0]

    if infmt == "sdf":
        mol = Chem.SDMolSupplier("../" + inlig, removeHs=False)[0]
    elif infmt == "mol2":
        mol = Chem.MolFromMol2File("../" + inlig, removeHs = False)
    runMethod(mol,basename,cut_type = "only1")
    os.chdir(datadir)
    
    ### run Vina ###
    datadir_frag = os.path.join(datadir, "Frags")
    if decoy:
        decoy_num = [line for line in os.path.join(datadir, "NumMol2.log")][0].rstrip()
        ### For CASF bechmark 
        run_Vina_Fragment(fn, inlig_pdb, datadir, datadir_frag, min = False, min_RW = False, RW = False, decoy = decoy)
        out = open(os.path.join(datadir,"NumFrags.csv"),"w")
        out.write("pdb,num_frag\n")
        for i in range(decoy_num + 1):
            i = str(i)
            ### i == 0 is for C ###
            ### i in [1,decoy_num] is for decoys fragments Vina score ###
            num_frags = generate_data(fn,i,datadir_frag)
        out.write(fn + "," + str(num_frags) + "\n")
        out.close()
    if water:
        run_Vina_Fragment(fn, inlig_pdb,datadir, datadir_frag, min = True, min_RW = True, RW = False, decoy = decoy)
        out = open(os.path.join(datadir,"NumFrags.csv"),"w")
        out.write("pdb,num_frag\n")
        for i in range(3):
            i = str(i)
            ### i == 0 is for C; i == 1 for Co; i == 2 for Crwo ###
            num_frags = generate_data(fn,i,datadir_frag)
        out.write(fn + "," + str(num_frags) + "\n")
        out.close()

    else:
        run_Vina_Fragment(fn, inlig_pdb, datadir, datadir_frag, min = True, min_RW = False, RW = False, decoy = False)
        out = open(os.path.join(datadir,"NumFrags.csv"),"w")
        out.write("pdb,num_frag\n")
        for i in range(2):
            i = str(i)
            ### i == 0 is for C; i == 1 for Co ###
            num_frags = generate_data(fn,i,datadir_frag)
        out.write(fn + "," + str(num_frags) + "\n")
        out.close()
    os.chdir(olddir)

    return None 


def get_input(datadir, fn):
    '''
    Get input files based on pdbid
    
    '''

    olddir = os.getcwd()
    os.chdir(datadir)
    ### get ligand ###
    ### input should be provided as either of sdf file (best choice) or mol2 file ###
    inlig1 = fn + "_ligand.mol2"
    inlig2 = fn + "_ligand.sdf"
    inlig3 = fn + "_ligand.pdb"

    ### check ligand input file ###
    if inlig2 in os.listdir("."):
        inlig = inlig2
        mol = Chem.SDMolSupplier(inlig,removeHs = False)[0]
        if mol == None:
            if inlig1 in os.listdir("."):
                inlig = inlig1
                mol = Chem.MolFromMol2File(inlig,removeHs = False)
                if mol == None:
                    error_massage = "Error:input ligand should be checked"
                    sys.exit(error_massage)
                else:
                    inlig_rdkit = inlig1
            else:
                error_massage = "Error:input ligand(sdf) should be checked"
                sys.exit(error_message)
        else:
            inlig_rdkit = inlig2

    else:
        inlig = inlig1
        mol = Chem.MolFromMol2File(inlig,removeHs = False)
        if mol == None:
            sys.exit("Error:input ligand(mol2) should be checked\ntry sdf")
        else:
            inlig_rdkit = inlig1

    if inlig3 not in os.listdir("."):
        inlig = inlig_rdkit
        infmt = inlig.split(".")[-1]
        outlig = inlig3
        cmd = obabel + " -i" + infmt + " " +  inlig + " -opdb -O " + outlig
        os.system(cmd)    

    ### get protein ###
    ### at least one protein structure should be provided with all waters ###
    inpro1 = fn + "_protein.pdb"
    inpro2 = fn + "_protein_all.pdb"

    if inpro1 not in os.listdir("."):
        inpro = inpro2
        outpro = open(inpro1,"w")
        protein_lines = get_pdbinfo.pdbinfo(fn, file = inpro).getProteinWaters()[0]
        
        outpro.write("".join(protein_lines))
        outpro.close()

    ### check input structures ###
    if os.path.isfile(inlig_rdkit) and os.stat(inlig_rdkit).st_size != 0:
        print("Ligand for conformation stability:" + inlig_rdkit)
    else:
        sys.exit("Error: rdkit input")
    if os.path.isfile(inlig3) and os.stat(inlig3).st_size != 0:
        print("Ligand for Vina, SASA, BA, ION:" + inlig3)  
    else:
        sys.exit("Error: ligand input (pdb)")
    if os.path.isfile(inpro1) and os.stat(inpro1).st_size != 0:
        print("Protein without water molecules:" + inpro1)
    else:
        sys.exit("Error: protein input without water")
    if os.path.isfile(inpro2) and os.stat(inpro2).st_size != 0:
        print("Protein with water molecules:" + inpro2)
    else:
        sys.exit("Error: protein input with water")
    
    os.chdir(olddir)

    return inlig_rdkit, inlig3, inpro1, inpro2



def run_features(datadir,fn, rw = False, water = True, decoy = False):
    '''
    run features

    '''


    inlig_rdkit, inlig_pdb, inpro_pro, inpro_water = get_input(datadir,fn)
    print("Finish Input Preparation")

    ### receptor water ###
    print("Receptor Water:" + str(rw))
    if rw:
        get_Crw(fn,inpro_pro,inpro_water,datadir)
        print("Finish generate RW")
    elif inpro_pro.split(".")[0] + "_RW.pdb" in os.listdir(datadir):
        print("Finish generate RW")
    else:
        olddir = os.getcwd()
        os.chdir(datadir)
        cmd = "cp " + inpro_water + " " + inpro_pro.split(".")[0] + "_RW.pdb"
        os.system(cmd)
        os.chdir(olddir)
        print("Finish copy RW")

    ### get Co, Crwo ###
    for st in ["","RW"]:
        get_Co(datadir,fn, inlig_pdb, st)
    print("Finish Co, Crwo")

    ### update input structures ###
    inlig_C = inlig_pdb
    inlig_Co = fn + "_lig_min.pdb"
    inlig_Crwo = fn + "_lig_min_RW.pdb"

    inpro_C = inpro_pro
    inpro_Crw = fn + "_protein_RW.pdb"

    ################################

    if water:
        d_type = {"":[inpro_C, inlig_C],"_min": [inpro_C, inlig_Co],"_min_RW": [inpro_Crw, inlig_Crwo]}
        ### get bw ###
        out_total = open(os.path.join(datadir,"Feature_BW_min_RW.csv"),"w")
        out_total.write("pdb,Nbw,Epw,Elw\n")
        cal_BW(out_total,fn,inpro_pro,inlig_Crwo,inpro_Crw,datadir)
        out_total.close()
        print("Finish BW")

    else:
        d_type = {"":[inpro_C, inlig_C],"_min": [inpro_C, inlig_Co]}
    
    for i in d_type.keys():
        inpro = d_type[i][0]
        inlig = d_type[i][1]
        ### get Vina58 ###
        outfile = open(os.path.join(datadir,"Vina58" + i + ".csv"),"w")
        outfile.write('pdb,vina,' + ','.join(['vina' + str(n+1) for n in range(58)]) + "\n")
        featureVina(outfile, fn, inpro, inlig, datadir)
        outfile.close()
        print("Finish Vina")

        ### get sasa ###
        f_type = ["P","N","DA","D","A","AR","H","PL","HA","SA"]
        f_feature = ["P2." + i for i in f_type] + ["P2dl." + i for i in f_type] + ["P2dp." + i for i in f_type]

        out_SASA = open(os.path.join(datadir, "SASA" + i + ".csv"),"w")
        out_SASA.write("pdb," + ",".join(f_feature) + "\n")
        cal_SASA(out_SASA,fn,inlig,inpro,datadir)
        out_SASA.close()
        print("Finish SASA")

        ### get ion ###
        outfile = open(os.path.join(datadir,"Num_Ions" + i + ".csv"),"w")
        outfile.write("pdb,Ni\n")
        cal_Ni(outfile,fn, inpro, inlig, datadir)
        outfile.close()
        print("Finish Ion")

    ### get dERMSD ###
    outfile = open(os.path.join(datadir,"dE_RMSD.csv"),"w")
    outfile.write("pdb,dE_global,RMSD_global,number0,number1\n")
    feature_cal(outfile,fn, inlig_rdkit, datadir, calc_type = "GenConfs")
    outfile.close()

    ### run fragments ###
    run_fragments(fn, datadir, inlig_rdkit, inlig_pdb, water = water, decoy = decoy)

    
    ### combine data ###
    if decoy: 
        combine(datadir,"", decoy = decoy)
    else:
        for i in d_type.keys():
            
            combine(datadir,i)
    
    print("Finish Feature Calculation")

    return None


if __name__ == "__main__":
    datadir = "/Users/jianinglu1/Documents/GitHub/deltaVinaXGB_develop/Test"
    fn = "01"
    run_features(datadir,fn)
