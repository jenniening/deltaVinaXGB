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

def renumber(fmt, infile, outfile):
    '''
    rename atoms in file based on order 

    '''

    if fmt == "mol2":
        lines = [line for line in open(infile)]
        atom_index = lines.index("@<TRIPOS>ATOM\n")
        bond_index = lines.index("@<TRIPOS>BOND\n")
        atom_lines = lines[atom_index + 1: bond_index]
        atoms = set([atom.split()[5].split(".")[0] for atom in atom_lines])
        atom_dic ={key:1 for key in atoms}
        atom_lines_new = []
        for line in atom_lines:
            atom_key = line.split()[5].split(".")[0]
            atom_old = line.split()[1]
            atom_new = atom_key.upper() + str(atom_dic[atom_key])
            if len(atom_old) > len(atom_new):
                atom_new = atom_new + (len(atom_old) -len(atom_new)) * " "
            elif len(atom_old) < len(atom_new):
                atom_old = atom_old + (len(atom_new) -len(atom_old)) * " "

            newline = line.replace(atom_old,atom_new,1)
            atom_lines_new.append(newline)
            atom_dic[atom_key] += 1

        outfile = open(outfile,"w")
        outfile.write("".join(lines[0:atom_index + 1]))
        outfile.write("".join(atom_lines_new))
        outfile.write("".join(lines[bond_index:]))
        outfile.close()
    elif fmt == "pdb":
        lines = [line for line in open(infile)]
        atom_lines = get_pdbinfo.pdbinfo(file = infile).getAtoms()
        atom_index = lines.index(atom_lines[0]) 
        bond_index = lines.index(atom_lines[-1]) + 1
        mol = Chem.MolFromPDBFile(infile,removeHs=False)
        atom_list = [atom.GetSymbol() for atom in mol.GetAtoms()]
        atoms = set([atom.GetSymbol() for atom in mol.GetAtoms()])
        atom_dic ={key:1 for key in atoms}
        atom_lines_new = []
        for idx, line in enumerate(atom_lines):
            atom_key = atom_list[idx]
            atom_old = get_pdbinfo.atmn(line).strip()
            atom_new = atom_key.upper() + str(atom_dic[atom_key])
            if len(atom_old) > len(atom_new):
                atom_new = atom_new + (len(atom_old) -len(atom_new)) * " "
            elif len(atom_old) < len(atom_new):
                atom_old = atom_old + (len(atom_new) -len(atom_old)) * " "
            newline = line.replace(atom_old,atom_new,1)
            atom_lines_new.append(newline)
            atom_dic[atom_key] += 1

        outfile = open(outfile,"w")
        outfile.write("".join(lines[0:atom_index]))
        outfile.write("".join(atom_lines_new))
        outfile.write("".join(lines[bond_index:]))
        outfile.close()

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

    ### correcting atom name for sasa calculation ###
    if inlig3 not in os.listdir("."):
        inlig = inlig_rdkit
        infmt = inlig.split(".")[-1]
        if infmt == "mol2":
            outlig_num = inlig.split(".")[0] + "_rename.mol2"
            renumber(infmt,inlig,outlig_num)
            outlig = inlig3.split(".")[0] + "_rename.pdb"
            cmd = obabel + " -i" + infmt + " " +  outlig_num + " -opdb -O " + outlig
            os.system(cmd)
            os.system("rm " + outlig_num)
        elif infmt == "sdf":
            outlig = inlig3
            cmd = obabel + " -i" + infmt + " " +  inlig + " -opdb -O " + outlig
            os.system(cmd)
            outlig_num  = outlig.split(".")[0] + "_rename.pdb"
            renumber('pdb', outlig, outlig_num)
    else:
        infmt = inlig3.split(".")[-1]
        outlig_num = inlig.split(".")[0] + "_rename.pdb"
        renumber(infmt,inlig3, outlig_num)

    inlig3 = fn + "_ligand_rename.pdb"

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


def write_decoys(infile, fn):
    '''
    split decoys into seperate file

    '''

    lines = open(infile).readlines()
    n = 0
    flag = False
    for line in lines:
        if line[0] != "#":
            if line == "@<TRIPOS>MOLECULE\n":
                flag = True
                n += 1
                if n == 1:
                    newfile = open(fn + "_" + str(n) + "_decoy.mol2","w")
                else:
                    newfile.close()
                    newfile = open(fn + "_" + str(n) + "_decoy.mol2","w")
            if flag:
                newfile.write(line)
    num_file = open("NumDecoys.csv","w")
    num_file.write(str(n))
    num_file.close()

    return n

def rewrite_decoy(ref_lig, decoy_list, ref_fmt, method = "order"):
    '''
    rewrite decoy based on ref_lig atom type and bond connection

    method: can be order (atom orders should be same)
            or name (atom names should be same)
    '''

    if ref_fmt == "mol2":
        if method == "order":
            for decoy in decoy_list:
                lines = open(decoy).readlines()
                atom_idx = lines.index("@<TRIPOS>ATOM")

            


    elif ref_fmt == "sdf":

    
def get_input_decoy(datadir, datadir_decoy, fn):

    inlig1 = fn + "_ligand.mol2"
    inlig2 = fn + "_ligand.sdf"

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


    ref_ligand = inlig_rdkit
    ref_fmt = ref_ligand.split(".")[1]
    if ref_fmt == "mol2":
        ref_ligand = ref_ligand.split(".")[0] + "_rename.mol2"
    ref_proten = fn + "_protein.pdb"
    decoy_file = fn + "_decoy.mol2"
    os.system("mkdir " + fn)
    os.chdir(fn)
    infile = os.path.join(datadir_decoy,fn + "_decoys.mol2")
    num = write_decoys(infile, fn)
    decoy_list = [fn + "_decoy.mol2" for i in range(1,num + 1)]
    rewrite_decoy(ref_ligand, decoy_list, ref_fmt)
    for decoy in decoy_list:
        if decoy_fmt == "mol2":
            obabel convert to pdb
        else:
            obable convert to sdf 



    return ref_lig, ref_protein, decoy_rdkit_list, decoy_list




def run_features(datadir,fn, f_type = "rw", rewrite = True, decoy = "CASF-docking"):
    '''

    run features
    
    f_type: features calculation type, defaults to "RW"
            "rw" --> get receptor water based our criteria and calculate score with receptor water 
            "w"  --> calculate score with water molecules in protein_all.pdb
            "nw" --> only calculate scores with no water (C and Co)
            "n"  --> only calculate scores for C
    
    rewrite: whether to rewrite all features

    decoy: CASF decoys or not

    '''
    if decoy:
        f_type = "n"
    else:
        inlig_rdkit, inlig_pdb, inpro_pro, inpro_water = get_input(datadir,fn)
        print("Finish Input Preparation")

    ### receptor water ###
    if f_type == "rw":
        print("Receptor Water: recalculate")
        if rewrite:
            get_Crw(fn,inpro_pro,inpro_water,datadir)
            print("Finish generate RW")
        else:
            if inpro_pro.split(".")[0] + "_RW.pdb" not in os.listdir(datadir):
                get_Crw(fn,inpro_pro,inpro_water,datadir)
                print("Finish generate RW")
            else:
                print("Finish generate RW")
    elif f_type == "w":
        print("Receptor Water: use waters in " + fn + "_protein_all.pdb" )
        if rewrite:
            olddir = os.getcwd()
            os.chdir(datadir)
            cmd = "cp " + inpro_water + " " + inpro_pro.split(".")[0] + "_RW.pdb"
            os.system(cmd)
            os.chdir(olddir)
            print("Finish copy RW")
        else:
            if inpro_pro.split(".")[0] + "_RW.pdb" not in os.listdir(datadir):
                olddir = os.getcwd()
                os.chdir(datadir)
                cmd = "cp " + inpro_water + " " + inpro_pro.split(".")[0] + "_RW.pdb"
                os.system(cmd)
                os.chdir(olddir)
                print("Finish copy RW")
            else:
                print("Finish generate RW")
    else:
        prinet("No Receptor Water")


    ### get Co, Crwo ###
    if f_type == "rw" or f_type == "w":
        for st in ["","RW"]:
            get_Co(datadir,fn, inlig_pdb, st)
        print("Finish Co, Crwo")
    elif f_type == "nw":
        for st in [""]:
            get_Co(datadir,fn, inlig_pdb, st)
        print("Finish Co")
    else:
        print("No Co")


    ### update input structures ###
    inlig_C = inlig_pdb
    inlig_Co = fn + "_lig_min.pdb"
    inlig_Crwo = fn + "_lig_min_RW.pdb"

    inpro_C = inpro_pro
    inpro_Crw = fn + "_protein_RW.pdb"

    ################################

    if f_type == "rw" or f_type == "w":
        d_type = {"":[inpro_C, inlig_C],"_min": [inpro_C, inlig_Co],"_min_RW": [inpro_Crw, inlig_Crwo]}
        ### get bw ###
        out_total = open(os.path.join(datadir,"Feature_BW_min_RW.csv"),"w")
        out_total.write("pdb,Nbw,Epw,Elw\n")
        cal_BW(out_total,fn,inpro_pro,inlig_Crwo,inpro_Crw,datadir)
        out_total.close()
        print("Finish BW")
    elif f_type == "nw"
        d_type = {"":[inpro_C, inlig_C],"_min": [inpro_C, inlig_Co]}
    else:
        d_type = {"":[inpro_C, inlig_C]}
         
    
    for i in d_type.keys():
        if i == "":
            print("C")
        elif i == "_min":
            print("Co")
        else:
            print("Crwo")

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
