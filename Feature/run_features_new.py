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

def run_fragments(fn, datadir, inlig, inlig_pdb, opt = None, decoy = False, decoy_list = None, decoy_type = None, decoy_pro = None):
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
    if infmt == "sdf":
        mol = Chem.SDMolSupplier("../" + inlig, removeHs=False)[0]
    elif infmt == "mol2":
        mol = Chem.MolFromMol2File("../" + inlig, removeHs = False)
    runMethod(mol,inlig,cut_type = "only1")
    os.chdir(datadir)
    
    ### run Vina ###
    datadir_frag = os.path.join(datadir, "Frags")

    if decoy:
        ### For CASF bechmark ###
        run_Vina_Fragment(fn, inlig_pdb, datadir, datadir_frag, decoy = True, decoy_list = decoy_list, decoy_pro = decoy_pro)
        out = open(os.path.join(datadir,"NumFrags_decoys" + decoy_type + ".csv"),"w")
        out.write("pdb,idx,num_frag\n")
        out_core = open(os.path.join(datadir,"Vina_core_decoys" + decoy_type + ".csv"),"w")
        out_core.write("pdb,idx,vina," + ','.join(['vina' + str(n+1) for n in range(58)]) + "\n")
        out_side = open(os.path.join(datadir,"Vina_side_decoys" + decoy_type + ".csv"),"w")
        out_side.write("pdb,idx,vina," + ','.join(['vina' + str(n+1) for n in range(58)]) + "\n")
        for idx, decoy in enumerate(decoy_list):
            idx_decoy = decoy.split("_")[1]
            ### i == 0 is for C ###
            ### i in [1,len(decoy_list+1)] is for decoys fragments Vina score ###
            num_frags = generate_data(fn,str(idx + 1),datadir_frag)
            lines = open(os.path.join(datadir_frag,"Vina_core_" + str(idx + 1) + ".csv")).readlines()
            out_core.write(fn + "," + idx_decoy + "," + ",".join(lines[1].split(",")[1:]))
            lines = open(os.path.join(datadir_frag,"Vina_side_" + str(idx + 1) + ".csv")).readlines()
            out_side.write(fn + "," + idx_decoy + "," + ",".join(lines[1].split(",")[1:]))
            out.write(fn + "," + idx_decoy + "," + str(num_frags) + "\n")
        out_core.close()
        out_side.close()
        out.close()
    else:
        ### for ligands ###
        ### all options in run_Vina_Fragment default to False/None ###
        header = "pdb,vina," + ','.join(['vina' + str(n+1) for n in range(58)]) + "\n"
        if opt == "rbwo":
            ### i == 0 is for C; i == 1 for Co; i == 2 for Crwo, i == 3 for Cbwo ###
            run_Vina_Fragment(fn, inlig_pdb, datadir, datadir_frag, min = True, min_RW = True, min_BW = True)
            name_type = ["","_min","_min_RW","_min_BW"]
        elif opt == "rwo":
            ### i == 0 is for C; i == 1 for Crwo ###
            run_Vina_Fragment(fn, inlig_pdb,datadir, datadir_frag,min_RW = True)
            name_type = ["","_min_RW"]
        elif opt == "bwo":
            ### i == 0 is for C; i == 1 for Cbwo ###
            run_Vina_Fragment(fn, inlig_pdb,datadir, datadir_frag,min_BW = True)
            name_type = ["","_min_BW"]
        elif opt == "pwo":
            ### i == 0 is for C; i == 1 for Cpwo ###
            run_Vina_Fragment(fn, inlig_pdb,datadir, datadir_frag,min_PW = True)
            name_type = ["","_min_PW"]
        elif opt == "o":
            ### i == 0 is for C; i == 1 for Co ###
            run_Vina_Fragment(fn, inlig_pdb,datadir, datadir_frag,min = True)
            name_type = ["","_min"]
        else:
            ### i == 0 is for C ###
            run_Vina_Fragment(fn, inlig_pdb,datadir, datadir_frag)
            name_type = [""]

        for i in range(len(name_type)):
            num_frags = generate_data(fn,str(i),datadir_frag)
            out_core = open(os.path.join(datadir,"Vina_core" + name_type[i] + ".csv"),"w")
            out_core.write(header)
            out_side = open(os.path.join(datadir,"Vina_side" + name_type[i] + ".csv"),"w")
            out_side.write(header)
            lines = open(os.path.join(datadir_frag,"Vina_core_" + str(i) + ".csv")).readlines()
            out_core.write(fn  + "," + ",".join(lines[1].split(",")[1:]))
            lines = open(os.path.join(datadir_frag,"Vina_side_" + str(i) + ".csv")).readlines()
            out_side.write(fn  + "," + ",".join(lines[1].split(",")[1:]))
        ### num_frag are same for all types of structure ###
        out = open(os.path.join(datadir,"NumFrags.csv"),"w")
        out.write("pdb,num_frag\n")
        out.write(fn + "," + str(num_frags) + "\n")
        out.close()

    os.chdir(datadir)

    if decoy:
        os.system("rm -r Frags")
    
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
    if inlig1 in os.listdir("."):
        inlig = inlig1
        mol = Chem.MolFromMol2File(inlig,removeHs = False)
        if mol == None:
            if inlig2 in os.listdir("."):
                inlig = inlig2
                mol = Chem.SDMolSupplier(inlig,removeHs = False)[0]
                if mol == None:
                    error_message = "Error:input ligand should be checked"
                    sys.exit(error_message)
                else:
                    inlig_rdkit = inlig2
            else:
                error_message = "Error:input ligand(mol2) should be checked"
                sys.exit(error_message)
        else:
            inlig_rdkit = inlig1
    else:
        inlig = inlig2
        mol = Chem.SDMolSupplier(inlig,removeHs = False)[0]
        if mol == None:
            sys.exit("Error:input ligand(sdf) should be checked\ntry sdf")
        else:
            inlig_rdkit = inlig2

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
            #os.system("rm " + outlig_num)
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


def write_decoys(infile, fn, ligname = False):
    '''
    split decoys into seperate file

    '''
    
    lines = open(infile).readlines()
    n = 0
    flag = False
    for idx, line in enumerate(lines):
        if line[0] != "#":
            if line == "@<TRIPOS>MOLECULE\n":
                if ligname:
                    name = lines[idx+1][0:4]
                    if name == fn:
                        flag = True
                        n += 1
                else:
                    flag = True
                    n += 1
                if n == 1:
                    newfile = open(fn + "_" + str(n) + "_decoy.mol2","w")
                elif n > 1:
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
        previous_lines = open(ref_lig).readlines()
        pre_atom_idx = previous_lines.index("@<TRIPOS>ATOM\n")
        pre_bond_idx = previous_lines.index("@<TRIPOS>BOND\n")
        ref_atoms = previous_lines[pre_atom_idx +1: pre_bond_idx]
        
        if method == "order":
            for decoy in decoy_list:
                lines = open(decoy).readlines()
                atom_idx = lines.index("@<TRIPOS>ATOM\n")
                bond_idx = lines.index("@<TRIPOS>BOND\n")
                atoms = lines[atom_idx +1:bond_idx]
                if len(ref_atoms) != len(atoms):
                    print("Num of atoms are not same: " + decoy)
                    continue
                atoms_coord = [ "%10s%10s%10s"%(line.split()[2],line.split()[3],line.split()[4]) for line in atoms]
                
                new_lines = [ line[0:17] + atoms_coord[idx] + line[46:] for idx, line in enumerate(ref_atoms)]
                newdecoy = open(decoy.split(".")[0] + "_correct.mol2","w")
                newdecoy.write("".join(lines[0:atom_idx + 1]))
                newdecoy.write("".join(new_lines))
                newdecoy.write("".join(previous_lines[pre_bond_idx:]))
                newdecoy.close()
        elif method == "name":
            ### nams shoule be unique ###
            for decoy in decoy_list:
                lines = open(decoy).readlines()
                atom_idx = lines.index("@<TRIPOS>ATOM\n")
                bond_idx = lines.index("@<TRIPOS>BOND\n")
                atoms = lines[atom_idx +1:bond_idx]
                if len(ref_atoms) != len(atoms):
                    print("Num of atoms are not same: " + decoy)
                    continue
                atoms_coord = {line[7:17].strip():"%10s%10s%10s"%(line.split()[2],line.split()[3],line.split()[4]) for line in atoms}
                new_lines = [line[0:17] + atoms_coord[line[7:17].strip()] + line[46:] for line in ref_atoms]
                newdecoy = open(decoy.split(".")[0] + "_correct.mol2","w")
                newdecoy.write("".join(lines[0:atom_idx + 1]))
                newdecoy.write("".join(new_lines))
                newdecoy.write("".join(previous_lines[pre_bond_idx:]))
                newdecoy.close()

    elif ref_fmt == "sdf":
        previous_lines = open(ref_lig).readlines()
        atom_num = int(previous_lines[3].split()[0:3])
        ref_atoms = previous_lines[4: atom_num+4]
        if method == "order":
            for decoy in decoy_list:
                lines = open(decoy).readlines()
                atom_idx = lines.index("@<TRIPOS>ATOM\n")
                bond_idx = lines.index("@<TRIPOS>BOND\n")
                atoms = lines[atom_idx +1:bond_idx]
                if len(ref_atoms) != len(atoms):
                    print("Num of atoms are not same: " + decoy)
                    continue
                atoms_coord = [ "%10s%10s%10s"%(line.split()[2],line.split()[3],line.split()[4]) for line in atoms]
                new_lines = [ atoms_coord[idx] + line[30:] for idx, line in enumerate(ref_atoms)]
                newdecoy = open(decoy.split(".")[0] + "_correct.sdf","w")
                newdecoy.write(lines[1])
                newdecoy.write("".join(previous_lines[1:4]))
                newdecoy.write("".join(new_lines))
                newdecoy.write("".join(previous_lines[atom_num+4:]))
                newdecoy.close()
        else:
            sys.exit("Error: ref structure should has same atom order with decoys if ref_fmt is sdf ")
    elif ref_fmt == "pdb":
        previous_lines = open(ref_lig).readlines()
        pre_atom_idx = [idx for idx, line in enumerate(previous_lines) if (line[0:6] == "ATOM  ") or (line[0:6] == "HETATM")][0]
        pre_bond_idx = [idx for idx, line in enumerate(previous_lines) if (line[0:6] == "ATOM  ") or (line[0:6] == "HETATM")][-1]
        ref_atoms = previous_lines[pre_atom_idx: pre_bond_idx + 1]
        if method == "order":
            for decoy in decoy_list:
                lines = open(decoy).readlines()
                atom_idx = lines.index("@<TRIPOS>ATOM\n")
                bond_idx = lines.index("@<TRIPOS>BOND\n")
                atoms = lines[atom_idx +1:bond_idx]
                if len(ref_atoms) != len(atoms):
                    print("Num of atoms are not same: " + decoy)
                    continue
                atoms_coord = [ "%8.3f%8.3f%8.3f"%(float(line.split()[2]),float(line.split()[3]),float(line.split()[4])) for line in atoms]
                new_lines = [ line[0:30] + atoms_coord[idx] + line[54:] for idx, line in enumerate(ref_atoms)]
                newdecoy = open(decoy.split(".")[0] + "_correct.pdb","w")
                newdecoy.write("COMPND    " + lines[1])
                newdecoy.write("AUTHOR    GENERATED BY OPEN BABEL 2.4.1\n")
                newdecoy.write("".join(new_lines))
                newdecoy.write("".join(previous_lines[pre_bond_idx + 1:]))
                newdecoy.close()
        elif method == "name":
            ### nams shoule be unique ###
            for decoy in decoy_list:
                lines = open(decoy).readlines()
                atom_idx = lines.index("@<TRIPOS>ATOM\n")
                bond_idx = lines.index("@<TRIPOS>BOND\n")
                atoms = lines[atom_idx +1:bond_idx]
                if len(ref_atoms) != len(atoms):
                    print("Num of atoms are not same: " + decoy)
                    continue
                atoms_coord = {line[7:17].strip():"%8.3f%8.3f%8.3f"%(float(line.split()[2]),float(line.split()[3]),float(line.split()[4])) for line in atoms}
                new_lines = [line[0:30] + atoms_coord[line[12:16].strip()] + line[54:] for line in ref_atoms]
                newdecoy = open(decoy.split(".")[0] + "_correct.pdb","w")
                newdecoy.write("COMPND    " + lines[1])
                newdecoy.write("AUTHOR    GENERATED BY OPEN BABEL 2.4.1")
                newdecoy.write("".join(new_lines))
                newdecoy.write("".join(previous_lines[pre_bond_idx + 1:]))
                newdecoy.close()



    return None

    
def get_input_decoy(datadir, datadir_decoy, fn, ligname):

    inlig1 = fn + "_ligand.mol2"
    inlig2 = fn + "_ligand.sdf"
    inlig3 = fn + "_ligand_rename.pdb"

    ### check ligand input file ###
    olddir = os.getcwd()
    os.chdir(datadir)
    if inlig1 in os.listdir("."):
        inlig = inlig1
        mol = Chem.MolFromMol2File(inlig,removeHs = False)
        if mol == None:
            if inlig2 in os.listdir("."):
                inlig = inlig2
                mol = Chem.SDMolSupplier(inlig,removeHs = False)[0]
                if mol == None:
                    error_message = "Error:input ligand should be checked"
                    sys.exit(error_message)
                else:
                    inlig_rdkit = inlig2
            else:
                error_message = "Error:input ligand(mol2) should be checked"
                sys.exit(error_message)
        else:
            inlig_rdkit = inlig1
    else:
        inlig = inlig2
        mol = Chem.SDMolSupplier(inlig,removeHs = False)[0]
        if mol == None:
            sys.exit("Error:input ligand(sdf) should be checked\ntry sdf")
        else:
            inlig_rdkit = inlig2
    os.chdir(olddir)

    ### prepare decoy file ####
    ref_ligand_rdkit = inlig_rdkit
    ref_ligand_pdb = inlig3
    ### copy reference ligand ###
    cmd = "cp " + os.path.join(datadir,ref_ligand_rdkit) + " " + datadir_decoy
    os.system(cmd)
    cmd = "cp " + os.path.join(datadir,ref_ligand_pdb) + " " + datadir_decoy
    os.system(cmd)
    ref_fmt = ref_ligand_rdkit.split(".")[1]
    decoy_file = fn + "_decoys.mol2"
    os.chdir(datadir_decoy)
    if decoy_file in os.listdir("../"):
        ### This is for CASF, which puts all decoys into one file ###
        os.system("cp ../" + fn + "_decoys.mol2 .")
        ### ligname True for CASF-2016 screening, since for target-ligand there are so many other structures in decoy mol2 file ### 
        num = write_decoys(decoy_file, fn, ligname = ligname)
    else:
        ### This is for user generated decoy files, default is 20, which is the Vina docking files generated by us ###
        num = 20
    decoy_list = [fn + "_" + str(i) + "_decoy.mol2" for i in range(1,num + 1)]
    rewrite_decoy(ref_ligand_rdkit, decoy_list, ref_fmt, "order")
    rewrite_decoy(ref_ligand_pdb, decoy_list, "pdb", "order")
    
    decoy_rdkit_list = [fn + "_" + str(i) + "_decoy_correct." + ref_fmt for i in range(1,num + 1) if fn + "_" + str(i) + "_decoy_correct." + ref_fmt in os.listdir(".")]
    decoy_pdb_list = [fn + "_" + str(i) + "_decoy_correct.pdb" for i in range(1,num + 1) if fn + "_" + str(i) + "_decoy_correct.pdb" in os.listdir(".")]
    if len(decoy_list) != len(decoy_rdkit_list):
        print(str(len(decoy_list)-len(decoy_rdkit_list)) + " decoys have been removed")
        print(decoy_rdkit_list)
    assert len(decoy_rdkit_list) == len(decoy_pdb_list), "Decoy File Prepration Failed"


    return ref_ligand_rdkit, ref_ligand_pdb, decoy_rdkit_list, decoy_pdb_list



def prepare_rw_receptor(datadir, fn, inpro_pro, inpro_water, inlig, water_type, rewrite = False):
    '''
    prepare protein with water for different water type

    water_type: protein part water type, defaults to "rbw"
            "rbw" --> get protein part water based RW and BW
            "rw" --> get protein part water based RW
            "bw" --> get protein part water based BW
            "pw" --> get protein part water based PW (all waters in protein_all.pdb are considered as protein part water)
            "n"  --> no consideration of water molecules

    '''
    if water_type == "rbw":
        print("Protein Water: recalculate by both RW and BW")
        if rewrite:
            get_Crw(fn,inpro_pro,inpro_water,datadir)
            print("Finish generate RW")
            out_total = open("Feature_BW_initial.csv","w")
            cal_BW(out_total,fn,inpro_pro,inlig,inpro_water,datadir, Feature = False)
            out_total.close()
            print("Finish generate BW")
        else:
            if inpro_pro.split(".")[0] + "_RW.pdb" not in os.listdir(datadir):
                get_Crw(fn,inpro_pro,inpro_water,datadir)
            else:
                print("Use previous generated RW")
            if inpro_pro.split(".")[0] + "_BW.pdb" not in os.listdir(datadir):
                out_total = open("Feature_BW_initial.csv","w")
                cal_BW(out_total,fn,inpro_pro,inlig,inpro_water,datadir, Feature = False)
                out_total.close()
                print("Finish generate BW")
            else:
                print("Use previous RW and BW")

    elif water_type == "rw":
        print("Protein Water: recalculate by RW")
        if rewrite:
            get_Crw(fn,inpro_pro,inpro_water,datadir)
            print("Finish generate RW")
        else:
            if inpro_pro.split(".")[0] + "_RW.pdb" not in os.listdir(datadir):
                get_Crw(fn,inpro_pro,inpro_water,datadir)
                print("Finish generate RW")
            else:
                print("Use previous RW")
    elif water_type == "pw":
        if rewrite:
            olddir = os.getcwd()
            os.chdir(datadir)
            cmd = "cp " + inpro_water + " " + inpro_pro.split(".")[0] + "_PW.pdb"
            os.system(cmd)
            os.chdir(olddir)
            print("Finish generated PW")
        else:
            if inpro_pro.split(".")[0] + "_RW.pdb" not in os.listdir(datadir):
                olddir = os.getcwd()
                os.chdir(datadir)
                cmd = "cp " + inpro_water + " " + inpro_pro.split(".")[0] + "_PW.pdb"
                os.system(cmd)
                os.chdir(olddir)
                print("Finish generated PW")
            else:
                print("Use previous PW")

    elif water_type == "bw":
        print("Protein Water: recalculate by BW")
        if rewrite:
            out_total = open("Feature_BW_initial.csv","w")
            cal_BW(out_total,fn,inpro_pro,inlig,inpro_water,datadir, Feature = False)
            out_total.close()
            print("Finish generate BW")
        else:
            if inpro_pro.split(".")[0] + "_BW.pdb" not in os.listdir(datadir):
                out_total = open("Feature_BW_initial.csv","w")
                cal_BW(out_total,fn,inpro_pro,inlig,inpro_water,datadir, Feature = False)
                out_total.close()
                print("Finish generate BW")
            else:
                print("Use previous BW")

    return None


def prepare_opt(datadir, fn, inlig_pdb, opt_type, decoy = False, pro = None, rewrite = False):
    '''
    prepare AutoDock Vina optimized ligand for different optimization type

    opt_type: receptor water type, defaults to "rbwo"
            "rbwo" --> do optimization for protein without water, with RW and with BW
            "rwo" --> do optimization for protein with RW
            "bwo" --> do optimization for protein with BW
            "pwo" --> do optimization for protein with PW(all waters in protein_all.pdb are considered as protein part water)
            "o" --> do optimization for protein without water
            "n" --> no optimization

    '''

    if opt_type == "rbwo":
        for st in ["","_RW","_BW"]:
            if rewrite:
                get_Co(datadir,fn, inlig_pdb, st, decoy = decoy, pro = pro)
            else:
                if fn + "_lig_min" + st + ".pdb" not in os.listdir(datadir):
                    get_Co(datadir,fn, inlig_pdb, st, decoy = decoy, pro = pro)
                else:
                    print("Use pervious generated C" + st + "O")
        print("Finish Optimization")

    elif opt_type == "rwo":
        for st in ["_RW"]:
            if rewrite:
                get_Co(datadir,fn, inlig_pdb, st, decoy = decoy, pro = pro)
            else:
                if fn + "_lig_min" + st + ".pdb" not in os.listdir(datadir):
                    get_Co(datadir,fn, inlig_pdb, st, decoy = decoy, pro = pro)
                else:
                    print("Use pervious generated C" + st + "O")
        print("Finish Optimization")

    elif opt_type == "bwo":
        for st in ["_BW"]:
            if rewrite:
                get_Co(datadir,fn, inlig_pdb, st, decoy = decoy, pro = pro )
            else:
                if fn + "_lig_min" + st + ".pdb" not in os.listdir(datadir):
                    get_Co(datadir,fn, inlig_pdb, st, decoy = decoy, pro = pro )
                else:
                    print("Use pervious generated C" + st + "O")
        print("Finish Optimization")

    elif opt_type == "pwo":
        for st in ["_PW"]:
            if rewrite:
                get_Co(datadir,fn, inlig_pdb, st, decoy = decoy, pro = pro)
            else:
                if fn + "_lig_min" + st + ".pdb" not in os.listdir(datadir):
                    get_Co(datadir,fn, inlig_pdb, st, decoy = decoy, pro = pro)
                else:
                    print("Use pervious generated C" + st + "O")
        print("Finish Optimization")

    elif opt_type == "o":
        for st in [""]:
            if rewrite:
                get_Co(datadir,fn, inlig_pdb, st, decoy = decoy, pro = pro)
            else:
                if (inlig_pdb.split(".")[0] + "_min.pdb" not in os.listdir(datadir)):
                    get_Co(datadir,fn, inlig_pdb, st, decoy = decoy, pro = pro)
                else:
                    print("Use pervious generated C" + st + "O")         
        print("Finish Optimization")
    
    return None


def prepare_opt_decoy(datadir_pro, datadir_decoy, fn, pro, decoy_pdb_list, opt_type ):
    '''
    optimize decoy structure provided by CASF-2013/CASF-2016 screening set

    '''

    os.chdir(datadir_decoy)
    ### copy protein structrue in current directory ###
    if opt_type == "rbwo":
        for st in ["","_RW","_BW"]:
            cmd = "cp " + os.path.join(datadir_pro, pro + "_protein" + st + ".pdb") + " ."
            os.system(cmd)
    elif opt_type == "rwo":
        cmd = "cp " + os.path.join(datadir_pro, pro + "_protein_RW.pdb") + " ."
        os.system(cmd)
    elif opt_type == "bwo":
        ### bw should be regenerated by currect ligand pose ###
        cmd = "cp " + os.path.join(datadir_pro, pro + "_protein_BW.pdb") + " ."
        os.system(cmd)
    elif opt_type == "pwo":
        cmd = "cp " + os.path.join(datadir_pro, pro + "_protein_PW.pdb") + " ."
        os.system(cmd)
    else:
        cmd = "cp " + os.path.join(datadir_pro, pro + "_protein.pdb") + " ."
        os.system(cmd)
    for decoy in decoy_pdb_list:
        idx = decoy.split("_")[1]
        prepare_opt(datadir_decoy, fn, decoy, opt_type, decoy = True, pro = pro)
        print("Decoy" + idx)

    return None
        
def prepare_rec_decoy(datadir_decoy, datadir_pro, pro, water_type):
    '''
    copy protein part for decoy structures 

    '''

    ref_protein = pro + "_protein.pdb"
    ref_protein_1 = None
    ref_protein_2 = None
    ref_protein_3 = None

    if water_type == "rbw":
        ref_protein_1 = pro + "_protein_RW.pdb"
        ref_protein_2 = pro + "_protein_BW.pdb"
    elif water_type == "rw":
        ref_protein_1 = pro + "_protein_RW.pdb"
    elif water_type == "bw":
        ref_protein_1 = pro + "_protein_BW.pdb"
    elif water_type == "pw":
        ref_protein_1 = pro + "_protein_PW.pdb"

    ### copy ref protein ###
    if ref_protein not in os.listdir(datadir_decoy):
        cmd = "cp " + os.path.join(datadir_pro,ref_protein) + " " + datadir_decoy
        os.system(cmd)
    if ref_protein_1 and ref_protein_1 not in os.listdir(datadir_decoy):
        cmd = "cp " + os.path.join(datadir_pro,ref_protein_1) + " " + datadir_decoy
        os.system(cmd)
    if ref_protein_2 and ref_protein_2 not in os.listdir(datadir_decoy):
        cmd = "cp " + os.path.join(datadir_pro,ref_protein_2) + " " + datadir_decoy
        os.system(cmd)
    
    return None
    

###### ---------------------- run feature functions ---------------------#######

def run_Vina_features(datadir, d_type, fn, inpro, inlig, decoy = False, decoy_list = None):
    '''
    calculate Vina features

    '''

    if decoy:
        outfile_V58 = open(os.path.join(datadir,"Vina58_decoys" + d_type + ".csv"),"w")
        outfile_V58.write('pdb,idx,vina,' + ','.join(['vina' + str(n+1) for n in range(58)]) + "\n")
        for decoy in decoy_list:
            idx = decoy.split("_")[1]
            outfile = open(os.path.join(datadir,"Vina58" + idx + ".csv"),"w")
            featureVina(outfile, fn, inpro, decoy, datadir)
            outfile.close()
            lines = open(os.path.join(datadir,"Vina58" + idx + ".csv")).readlines()
            outfile_V58.write(fn + "," + idx + "," + ",".join(lines[0].split(",")[1:]) )
            rm_cmd = "rm " + os.path.join(datadir,"Vina58" + idx + ".csv")
            os.system(rm_cmd)
            print("Finish Vina" + idx)

    else:
        outfile = open(os.path.join(datadir,"Vina58" + d_type + ".csv"),"w")
        outfile.write('pdb,vina,' + ','.join(['vina' + str(n+1) for n in range(58)]) + "\n")
        featureVina(outfile, fn, inpro, inlig, datadir)
        outfile.close()
        print("Finish Vina")

    return None


def run_SASA_features(datadir, d_type, fn, inpro, inlig, decoy = False, decoy_list = None):
    '''
    calculate SASA features

    '''

    f_type = ["P","N","DA","D","A","AR","H","PL","HA","SA"]
    f_feature = ["P2." + i for i in f_type] + ["P2dl." + i for i in f_type] + ["P2dp." + i for i in f_type]
    if decoy:
        outfile_SASA = open(os.path.join(datadir,"SASA_decoys" + d_type + ".csv"),"w")
        outfile_SASA.write("pdb,idx," + ",".join(f_feature) + "\n")
        for decoy in decoy_list:
            idx = decoy.split("_")[1]
            outfile = open(os.path.join(datadir, "SASA" + idx + ".csv"),"w")
            try:
                cal_SASA(outfile,fn,decoy,inpro,datadir)
            except:
                ### if SASA failed, means there is no overlap between protein and decoy ###
                outfile.write(fn + "," + ",".join([str(0.00) for i in range(30)]) + "\n")
            outfile.close()
            lines = open(os.path.join(datadir, "SASA" + idx + ".csv")).readlines()
            outfile_SASA.write(fn + "," + idx + "," + ",".join(lines[0].split(",")[1:]))
            rm_cmd = "rm " + os.path.join(datadir,"SASA" + idx + ".csv")
            os.system(rm_cmd)
            print("Finish SASA" + idx)
    else:
        out_SASA = open(os.path.join(datadir, "SASA" + d_type + ".csv"),"w")
        out_SASA.write("pdb," + ",".join(f_feature) + "\n")
        cal_SASA(out_SASA,fn,inlig,inpro,datadir)
        out_SASA.close()
        print("Finish SASA")

    return None

def run_Ion_features(datadir, d_type, fn, inpro, inlig, decoy = False, decoy_list = None):
    '''
    calculate Ion features

    '''

    if decoy:
        outfile_Ion = open(os.path.join(datadir,"Num_Ions_decoys" + d_type + ".csv"),"w")
        outfile_Ion.write("pdb,idx,Ni\n")
        for decoy in decoy_list:
            idx = decoy.split("_")[1]
            outfile = open(os.path.join(datadir,"Num_Ions" + idx + ".csv"),"w")
            cal_Ni(outfile, fn, inpro, decoy, datadir)
            outfile.close()
            lines = open(os.path.join(datadir,"Num_Ions" + idx + ".csv")).readlines()
            outfile_Ion.write(fn + "," + idx  + "," + ",".join(lines[0].split(",")[1:]))
            rm_cmd = "rm " + os.path.join(datadir,"Num_Ions" + idx + ".csv")
            os.system(rm_cmd)
            print("Finish Ion" + idx)

    else:
        outfile = open(os.path.join(datadir,"Num_Ions" + d_type + ".csv"),"w")
        outfile.write("pdb,Ni\n")
        cal_Ni(outfile,fn, inpro, inlig, datadir)
        outfile.close()
        print("Finish Ion")

    return None


def run_BW_features(datadir, d_type, fn, inpro_pro, inlig, inpro, decoy = False, decoy_list = None):
    '''
    calculate BW features 

    '''

    if decoy:
        outfile_BW = open(os.path.join(datadir,"Feature_BW_decoys" + d_type + ".csv"),"w")
        outfile_BW.write("pdb,idx,Nbw,Epw,Elw\n")
        for decoy in decoy_list:
            idx = decoy.split("_")[1]
            idx = decoy.split("_")[1]
            outfile = open(os.path.join(datadir,"Feature_BW" + idx + ".csv"),"w")
            cal_BW(outfile,fn,inpro_pro,decoy,inpro,datadir)
            outfile.close()
            lines = open(os.path.join(datadir,"Feature_BW" + idx + ".csv")).readlines()
            outfile_BW.write(fn + "," + idx + "," + ",".join(lines[0].split(",")[1:]))
            rm_cmd = "rm " + os.path.join(datadir,"Feature_BW" + idx + ".csv")
            os.system(rm_cmd)
            print("Finish BW" + idx )
        outfile_BW.close()
    else:
        out_total = open(os.path.join(datadir,"Feature_BW" + d_type + ".csv"),"w")
        out_total.write("pdb,Nbw,Epw,Elw\n")
        cal_BW(out_total,fn,inpro_pro,inlig,inpro,datadir)
        out_total.close()
        print("Finish Bridging Water Feature Calculation")

    return None

def run_dE_features(datadir, fn, inlig_rdkit, decoy = False, decoy_rdkit_list = None, rewrite = False):
    '''
    calculate dE features

    '''

    if decoy:
        outfile_dE = open(os.path.join(datadir,"dE_RMSD_decoys.csv"),"w")
        outfile_dE.write("pdb,idx,dE_global,RMSD_global,number0,number1\n")
        for decoy in decoy_rdkit_list:
            idx = decoy.split("_")[1]
            outfile = open(os.path.join(datadir, "dE_RMSD" + idx + ".csv"),"w")
            feature_cal(outfile, fn, decoy, datadir, calc_type = "none", rewrite = rewrite)
            outfile.close()
            lines = open(os.path.join(datadir, "dE_RMSD" + idx + ".csv")).readlines()
            outfile_dE.write(fn + "," + idx + "," + ",".join(lines[0].split(",")[1:]))
            rm_cmd = "rm " + os.path.join(datadir,"dE_RMSD" + idx + ".csv")
            os.system(rm_cmd)
            print("Finish dE_RMSD" + idx)
        outfile_dE.close()
    else:
        outfile = open(os.path.join(datadir,"dE_RMSD.csv"),"w")
        outfile.write("pdb,dE_global,RMSD_global,number0,number1\n")
        feature_cal(outfile,fn, inlig_rdkit, datadir, calc_type = "GenConfs", rewrite = rewrite)
        outfile.close()

    return None

def run_Frag_features(datadir, fn, inlig_rdkit, inlig_pdb, opt_type, decoy = False, decoy_list = None, decoy_type = None, decoy_pro = None):
    '''
    calculate fragments features

    '''

    if decoy:
        run_fragments(fn, datadir, inlig_rdkit, inlig_pdb, decoy = True, decoy_list = decoy_list, decoy_type = decoy_type, decoy_pro = decoy_pro)
    else:
        run_fragments(fn, datadir, inlig_rdkit, inlig_pdb, opt = opt_type)

    return None


###### ---------------------- run feature calculations ---------------------#######    
def feature_calculation_decoy(datadir, datadir_pro, datadir_decoy, fn, pro, ref_ligand_rdkit,ref_ligand_pdb,decoy_rdkit_list,decoy_pdb_list,water_type = "n",opt_type = "n", feature_type = "all"):
    '''
    feature calculation for decoys

    for docking decoys: only has water_type but not opt_type
    for screeening decoys: has both water_type and opt_type, they are corresponding to each other

    '''
    ### update input structures ###
    inlig_C = decoy_pdb_list ### initial structure ###
    inlig_Co = [fn + "_" + decoy.split("_")[1] + "_decoy_min.pdb" for decoy in decoy_pdb_list]
    inlig_Crwo = [fn + "_" + decoy.split("_")[1] + "_decoy_min_RW.pdb" for decoy in decoy_pdb_list]
    inlig_Cbwo = [fn + "_" + decoy.split("_")[1] + "_decoy_min_BW.pdb" for decoy in decoy_pdb_list]
    inlig_Cpwo = [fn + "_" + decoy.split("_")[1] + "_decoy_min_PW.pdb" for decoy in decoy_pdb_list]

    inpro_C = pro + "_protein.pdb" ### protein without water ####
    inpro_Crw = pro + "_protein_RW.pdb"
    inpro_Cbw = pro + "_protein_BW.pdb"
    inpro_Cpw = pro + "_protein_PW.pdb"

    ################################
    if opt_type == "rbwo":
        d_type = {"":[inpro_C, inlig_C],"_min": [inpro_C, inlig_Co],"_min_RW": [inpro_Crw, inlig_Crwo], "_min_BW":[inpro_Cbw, inlig_Cbwo]}
    elif opt_type == "rwo":
        d_type = {"_min_RW": [inpro_Crw, inlig_Crwo]}
    elif opt_type == "bwo":
        d_type = {"_min_BW": [inpro_Cbw, inlig_Cbwo]}
    elif opt_type == "pwo":
        d_type = {"_min_PW":[inpro_Cpw, inlig_Cpwo]}
    elif opt_type == "o":
        d_type = {"_min":[inpro_C, inlig_Co]}
    else:
        if water_type == "rbw":
            d_type = {"":[inpro_C, inlig_C], "_RW":[inpro_Crw, inlig_C], "_BW":[inpro_Cbw, inlig_C]}
        elif water_type == "rw":
            d_type = {"_RW":[inpro_Crw, inlig_C]}
        elif water_type == "bw":
            d_type = {"_BW":[inpro_Cbw, inlig_C]}
        elif water_type == "pw":
            d_type = {"_BW":[inpro_Cpw, inlig_C]}
        else:
            d_type = {"":[inpro_C, inlig_C]}


    for i in d_type.keys():
        if i == "":
            print("C")
        elif i == "_min":
            print("Co")
        elif i == "_min_RW":
            print("Crwo")
        elif i == "_RW":
            print("Crw")
        elif i == "_min_BW":
            print("Cbwo")
        elif i == "_BW":
            print("Cbw")
        elif i == "_min_PW":
            print("Cpwo")
        elif i == "_PRW":
            print("Cpw")
        
        inpro = d_type[i][0]
        decoy_list = d_type[i][1]
        inlig = None
        inlig_rdkit = None
        ### get BW features ###
        ### only for structure with water in protein ###
        if i in ["_RW","_min_RW","_BW","_min_BW"]:
            inpro_pro = inpro_C
            if feature_type == "all" or feature_type == "BW":
                run_BW_features(datadir_decoy, i, fn, inpro_pro, inlig, inpro, decoy = True, decoy_list = decoy_list)
                print("Finish Bridging Water Feature Calculation")
            else:
                print("Use previous calculated Bridging Water Feature")

        ### get Vina58 ###
        if feature_type == "all" or feature_type == "Vina":
            run_Vina_features(datadir_decoy, i, fn, inpro, inlig, decoy = True, decoy_list = decoy_list)
            print("Finish Vina")
        else:
            print("Use previous calculated Vina")

        ### get sasa ###
        if feature_type == "all" or feature_type == "SASA":
            run_SASA_features(datadir_decoy, i, fn, inpro, inlig, decoy = True, decoy_list = decoy_list)
            print("Finish SASA")
        else:
            print("Use previous calculated SASA")

        ### get ion ###
        if feature_type == "all" or feature_type == "Ion":
            run_Ion_features(datadir_decoy, i, fn, inpro, inlig, decoy = True, decoy_list = decoy_list)
            print("Finish Ion")
        else:
            print("Use previous calculated Ion")

        ### run fragments ###
        if feature_type == "all" or feature_type == "Frag":
            run_Frag_features(datadir_decoy, fn, ref_ligand_rdkit, ref_ligand_pdb, opt_type = None, decoy = True, decoy_list = decoy_list, decoy_type = i, decoy_pro = inpro)

        else:
            print("Use previous calculated Frags")


    ### get dERMSD ###
    if feature_type == "all" or feature_type == "dE":
        if "dE_RMSD_decoys.csv" not in os.listdir(datadir_decoy):
            ### copy previous generated confs by other type decoys ###
            confs = os.path.join(datadir, fn + "_ligand_confs.sdf")
            lowest = os.path.join(datadir, fn + "_ligand_global_min.sdf")
            cmd = "cp " + confs + " " + datadir_decoy
            os.system(cmd)
            cmd = "cp " + lowest + " " + datadir_decoy
            run_dE_features(datadir_decoy, fn, inlig_rdkit, decoy = True, decoy_rdkit_list = decoy_rdkit_list)
        else:
            print("Use previous calculated dE")
    else:
        print("Use previous calculated dE")

    ### combine data ###
    for i in d_type.keys():
        combine(datadir_decoy,i,decoy = True)

    return None


def feature_calculation_ligand(datadir,fn, inlig_pdb, inlig_rdkit, inpro_pro, water_type, opt_type, rewrite = False, feature_type = "all"):
    '''
    feature calculation for ligands (C, Co, Crwo, Cbwo, Cpwo)
    
    '''

    ### update input structures ###
    inlig_C = inlig_pdb ### initial structure ###
    inlig_Co = fn + "_lig_min.pdb"
    inlig_Crwo = fn + "_lig_min_RW.pdb"
    inlig_Cbwo = fn + "_lig_min_BW.pdb"
    inlig_Cpwo = fn + "_lig_min_PW.pdb"

    inpro_C = inpro_pro ### protein without water ####
    inpro_Crw = fn + "_protein_RW.pdb"
    inpro_Cbw = fn + "_protein_BW.pdb"
    inpro_Cpw = fn + "_protein_PW.pdb"

    ################################

    ### get output file type ###
    if opt_type == "rbwo":
        d_type = {"":[inpro_C, inlig_C],"_min": [inpro_C, inlig_Co],"_min_RW": [inpro_Crw, inlig_Crwo], "_min_BW":[inpro_Cbw, inlig_Cbwo]}
    elif opt_type == "rwo":
        d_type = {"_min_RW": [inpro_Crw, inlig_Crwo]}
    elif opt_type == "bwo":
        d_type = {"_min_BW": [inpro_Cbw, inlig_Cbwo]}
    elif opt_type == "pwo":
        d_type = {"_min_PW":[inpro_Cpw, inlig_Cpwo]}
    elif opt_type == "o":
        d_type = {"_min":[inpro_C, inlig_Co]}
    else:
        d_type = {"":[inpro_C, inlig_C]}
         
    
    for i in d_type.keys():
        if i == "":
            print("C")
        elif i == "_min":
            print("Co")
        elif i == "_min_RW":
            print("Crwo")
        elif i == "_min_BW":
            print("Cbwo")
        elif i == "_min_PW":
            print("Cpwo")

        inpro = d_type[i][0]
        inlig = d_type[i][1]

        ### get BW features ###
        ### only for structure with water in protein ###

        if i in ["_min_RW", "_min_BW", "_min_PW"]:
            inpro_pro = inpro_C
            if feature_type == "all" or feature_type == "BW":
                run_BW_features(datadir, i, fn, inpro_pro, inlig, inpro)
                print("Finish Bridging Water Feature Calculation")
            else:
                print("Use previous calculated Bridging Water Feature")

        ### get Vina58 ###
        if feature_type == "all" or feature_type == "Vina":
            run_Vina_features(datadir, i, fn, inpro, inlig)
            print("Finish Vina")
        else:
            print("Use previous calculated Vina")

        ### get sasa ###
        if feature_type == "all" or feature_type == "SASA":
            run_SASA_features(datadir, i, fn, inpro, inlig)
            print("Finish SASA")
        else:
            print("Use previous calculated SASA")

        ### get ion ###
        if feature_type == "all" or feature_type == "Ion":
            run_Ion_features(datadir, i, fn, inpro, inlig)
            print("Finish Ion")
        else:
            print("Use previous calculated Ion")

    ### get dERMSD ###
    if feature_type == "all" or feature_type == "dE":
        run_dE_features(datadir, fn, inlig_rdkit, rewrite)
    else:
        print("Use previous calculated dE")
        
    ### run fragments ###
    if feature_type == "all" or feature_type == "Frag":
        run_Frag_features(datadir, fn, inlig_rdkit, inlig_pdb, opt_type)
    else:
        print("Use previous calculated Frags")

    ### combine data ###
    for i in d_type.keys():
        combine(datadir,i)
    
    print("Finish Feature Calculation")

    return None


def run_features(datadir, datadir_pro, datadir_decoy, fn, pro, water_type = "rbw", opt_type = "rbwo", decoy_type = "docking", rewrite = False, decoy = False, ligname = False, feature_type = "all"):
    '''

    run features
    
    datadir: directory for structures
    datadir_pro: directory for protein structure, only for decoys in screening test
    datadir_decoy: directory for decoy structures, only for decoys
    fn: pdbid
    pro: pdbid for protein, only for decoys in screening test

    water_type: protein part water type, defaults to "rbw"
            "rbw" --> get protein part water based RW and BW
            "rw" --> get protein part water based RW
            "bw" --> get protein part water based BW
            "pw" --> get protein part water based PW (all waters in protein_all.pdb are considered as protein part water)
            "n"  --> no consideration of water molecules

    opt_type: receptor water type, defaults to "rbwo"
            "rbwo" --> do optimization for protein without water, with RW and with BW
            "rwo" --> do optimization for protein with RW
            "bwo" --> do optimization for protein with BW
            "pwo" --> do optimization for protein with PW(all waters in protein_all.pdb are considered as protein part water)
            "o" --> do optimization for protein without water
            "n" --> no optimization

            
    rewrite: whether to rewrite all features, defaults to False

    decoy: CASF decoys or not, defaults to False

    ligname: whether use ligname to grab decoys (Only in CASF-2016_Screening), defaults to False

    feature_type: which feature will be calcuated, defaults to "all"
            "Vina" --> only calculate Vina58 features
            "SASA" --> only calculate SASA features
            "BW" --> only calculate BW features
            "dE" --> only calculate ligand satbility features
            "Ion" --> only calculate Ion features
            "Frag" --> only calculate Fragments features

    '''

    if decoy:
        if datadir_pro == None:
            datadir_pro = datadir
        if pro == None:
            pro = fn
        ### CASF-2013/2016 docking/screening, no water has been used in decoy preparation ###
        ### for CASF-2013/2016 docking, we can do also do local optimization for decoys with or without water molecules. And we need to recalculate the RMSD after finish ###
        ### previously, I think CASF-2013/2016 decoys scores should without water, but might be added water effect ###
        #water_type = "n"
        ref_ligand_rdkit, ref_ligand_pdb, decoy_rdkit_list, decoy_pdb_list = get_input_decoy(datadir,datadir_decoy, fn, ligname)
            
    else:
        inlig_rdkit, inlig_pdb, inpro_pro, inpro_water = get_input(datadir,fn)
        print("Finish Input Preparation")

    ### receptor water ###
    if water_type != "n":
        if not decoy:
            prepare_rw_receptor(datadir, fn, inpro_pro, inpro_water, inlig_pdb, water_type, rewrite = rewrite)
        else:
            prepare_rec_decoy(datadir_decoy, datadir_pro, pro, water_type)
            print("Consideration of Water Effect for Decoys")
    
    else:
        print("No Consideration of Water")

    ### get Co, Crwo ###
    if opt_type != "n":
        if not decoy:
            prepare_opt(datadir, fn, inlig_pdb, opt_type, rewrite = rewrite)
        else:
            print("Do optimization for decoys")
            prepare_opt_decoy(datadir_pro, datadir_decoy, fn, pro, decoy_pdb_list, opt_type)
    else:
        print("No Optimized Ligand")


    ### update input structures ###
    if decoy:
        feature_calculation_decoy(datadir,datadir_pro,datadir_decoy, fn, pro, ref_ligand_rdkit,ref_ligand_pdb,decoy_rdkit_list,decoy_pdb_list,
                                    water_type = water_type, opt_type = opt_type, feature_type = feature_type)
    else:
        feature_calculation_ligand(datadir,fn, inlig_pdb, inlig_rdkit, inpro_pro, water_type, opt_type, rewrite, feature_type)
    
    print("Finish Feature Calculation")

    return None


if __name__ == "__main__":
    datadir = "/Users/jianinglu1/Documents/GitHub/deltaVinaXGB_develop/Test"
    fn = "01"
    
 
