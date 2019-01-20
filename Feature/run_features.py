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

def run_fragments(fn, datadir, inlig, inlig_pdb, opt = True, water = True, decoy = False, decoy_list = None, decoy_type = None, decoy_pro = None):
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
        ### For CASF bechmark 
        run_Vina_Fragment(fn, inlig_pdb, datadir, datadir_frag, min = False, min_RW = False, RW = False, decoy = True, decoy_list = decoy_list, decoy_pro = decoy_pro)
        out = open(os.path.join(datadir,"NumFrags" + decoy_type + ".csv"),"w")
        out.write("pdb,idx, num_frag\n")
        out_core = open("Vina_core" + decoy_type + ".csv","w")
        out_side = open("Vina_side" + decoy_type + ".csv","w")
        for idx, decoy in enumerate(decoy_list):
            ### i == 0 is for C ###
            ### i in [1,len(decoy_list+1)] is for decoys fragments Vina score ###
            num_frags = generate_data(fn,idx + 1,datadir_frag)
            lines = open("Vina_core" + str(idx + 1) + "csv").readlines()
            out_core.write(fn + "," + str(idx +1) + "," + ",".join(lines[0].split(",")[1:]))
            lines = open("Vina_side" + str(idx + 1) + "csv").readlines()
            out_side.write(fn + "," + str(idx +1) + "," + ",".join(lines[0].split(",")[1:]))
            out.write(fn + "," + str(idx+1) + "," + str(num_frags) + "\n")
        out_core.close()
        out_side.close()
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

    elif water == False and opt == True
        run_Vina_Fragment(fn, inlig_pdb, datadir, datadir_frag, min = True, min_RW = False, RW = False, decoy = False)
        out = open(os.path.join(datadir,"NumFrags.csv"),"w")
        out.write("pdb,num_frag\n")
        for i in range(2):
            i = str(i)
            ### i == 0 is for C; i == 1 for Co ###
            num_frags = generate_data(fn,i,datadir_frag)
        out.write(fn + "," + str(num_frags) + "\n")
        out.close()
    elif water == False and opt == False:
        run_Vina_Fragment(fn, inlig_pdb, datadir, datadir_frag, min = False, min_RW = False, RW = False, decoy = False)
        out = open(os.path.join(datadir,"NumFrags.csv"),"w")
        out.write("pdb,num_frag\n")
        for i in range(1):
            i = str(i)
            ### i == 0 is for C ###
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

            newline = line.replace(atom_old,atom_new)
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
            newline = line.replace(atom_old,atom_new)
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
                    error_message = "Error:input ligand should be checked"
                    sys.exit(error_message)
                else:
                    inlig_rdkit = inlig1
            else:
                error_message = "Error:input ligand(sdf) should be checked"
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
                atoms_coord = {line[7:17].strip():"%10s%10s%10s"%(line.split()[2],line.split()[3],line.split()[4]) for line in atoms}
                new_lines = [line[0:17] + atoms_coord[line[7:17].strip()] + line[46:] for line in ref_atoms]
                newdecoy = open(decoy.split(".")[0] + "_correct.mol2","w")
                newdecoy.write("".join(lines[0:atom_idx + 1]))
                newdecoy.write("".join(new_lines))
                newdecoy.write("".join(previous_lines[pre_bond_idx:]))
                newdecoy.close()

    elif ref_fmt == "sdf":
        previous_lines = open(ref_lig).readlines()
        atom_num = previous_lines[3].split()[0]
        ref_atoms = previous_lines[4: atom_num+4]
        if method == "order":
            for decoy in decoy_list:
                lines = open(decoy).readlines()
                atom_idx = lines.index("@<TRIPOS>ATOM\n")
                bond_idx = lines.index("@<TRIPOS>BOND\n")
                atoms = lines[atom_idx +1:bond_idx]
                atoms_coord = [ "%10s%10s%10s"%(line.split()[2],line.split()[3],line.split()[4]) for line in atoms]
                new_lines = [ atoms_coord[idx] + line[30:] for idx, line in enumerate(ref_atoms)]
                newdecoy = open(decoy.split(".")[0] + "_correct.sdf","w")
                newdecoy.write("".join(lines[0:4]))
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
                atoms_coord = [ "%8s%8s%8s"%(line.split()[2][0:-1],line.split()[3][0:-1],line.split()[4][0:-1]) for line in atoms]
                new_lines = [ line[0:31] + atoms_coord[idx] + line[55:] for idx, line in enumerate(ref_atoms)]
                newdecoy = open(decoy.split(".")[0] + "_correct.pdb","w")
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
                atoms_coord = {line[7:17].strip():"%8s%8s%8s"%(line.split()[2][0:-1],line.split()[3][0:-1],line.split()[4][0:-1]) for line in atoms}
                new_lines = [line[0:16] + atoms_coord[line[12:16].strip()] + line[54:] for line in ref_atoms]
                newdecoy = open(decoy.split(".")[0] + "_correct.pdb","w")
                newdecoy.write("".join(lines[0:atom_idx + 1]))
                newdecoy.write("".join(new_lines))
                newdecoy.write("".join(previous_lines[pre_bond_idx:]))
                newdecoy.close()



    return None

    
def get_input_decoy(datadir, datadir_decoy, fn):

    inlig1 = fn + "_ligand.mol2"
    inlig2 = fn + "_ligand.sdf"
    inlig3 = fn + "_ligand_rename.pdb"

    ### check ligand input file ###
    if inlig2 in os.listdir("."):
        inlig = inlig2
        mol = Chem.SDMolSupplier(inlig,removeHs = False)[0]
        if mol == None:
            if inlig1 in os.listdir("."):
                inlig = inlig1
                mol = Chem.MolFromMol2File(inlig,removeHs = False)
                if mol == None:
                    error_message = "Error:input ligand should be checked"
                    sys.exit(error_message)
                else:
                    inlig_rdkit = inlig1
            else:
                error_message = "Error:input ligand(sdf) should be checked"
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

    ### prepare decoy file ####
    ref_ligand_rdkit = inlig_rdkit
    ref_ligand_pdb = inlig3
    ref_fmt = ref_ligand_rdkit.split(".")[1]
    decoy_file = fn + "_decoy.mol2"
    os.system("mkdir " + fn)
    os.chdir(fn)
    os.system("cp ../" + fn + "_decoy.mol2 .")
    num = write_decoys(decoy_file, fn)
    decoy_list = [fn + "_" + str(i) + "_decoy.mol2" for i in range(1,num + 1)]
    rewrite_decoy(ref_ligand_rdkit, decoy_list, ref_fmt, "order")
    rewrite_decoy(ref_ligand_pdb, decoy_list, "pdb", "order")
    os.file.st_size()
    decoy_rdkit_list = [fn + "_" + str(i) + "_decoy_correct." + ref_fmt for i in range(1,num + 1) if os.stat(fn + "_" + str(i) + "_decoy_correct." + ref_fmt).st_size != 0]
    decoy_list = [fn + "_" + str(i) + "_decoy_correct.pdb" for i in range(1,num + 1) if os.stat(fn + "_" + str(i) + "_decoy_correct.pdb").st_size != 0]

    assert len(decoy_rdkit_list) == len(decoy_list), "Decoy File Prepration Failed"


    return ref_ligand_rdkit, ref_ligand_pdb, decoy_rdkit_list, decoy_list



def prepare_rw_receptor(datadir, fn, inpro_pro, inpro_water, water_type, rewrite):
    '''
    prepare receptor water file for different water type

    water_type: receptor water type, defaults to "rw"
            "rw" --> get receptor water based our criteria
            "w"  --> all waters in protein_all.pdb are considered as receptor water
            "n"  --> no consideration of water molecules

    '''


    if water_type == "rw":
        print("Receptor Water: recalculate")
        if rewrite:
            get_Crw(fn,inpro_pro,inpro_water,datadir)
            print("Finish generate RW")
        else:
            if inpro_pro.split(".")[0] + "_RW.pdb" not in os.listdir(datadir):
                get_Crw(fn,inpro_pro,inpro_water,datadir)
                print("Finish generate RW")
            else:
                print("Use previous generated RW")
    elif water_type == "w":
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
                print("Use previous copy RW")

    return None


def prepare_opt_ligand(datadir, fn, inlig_pdb, opt_type, rewrite):
    if opt_type == "wo":
        for st in ["","RW"]:
            if rewrite:
                get_Co(datadir,fn, inlig_pdb, st)
            else:
                if (inlig_pdb.split(".")[0] + "_min.pdb" not in os.listdir(datadir) ) or (inlig_pdb.split(".")[0] + "_min_RW.pdb" not in os.listdir(datadir) ):
                    get_Co(datadir,fn, inlig_pdb, st)
                else:
                    print("Use pervious generated Co, Crwo")
        print("Finish Co, Crwo")

    elif opt_type == "o":
        for st in [""]:
            if rewrite:
                get_Co(datadir,fn, inlig_pdb, st)
            else:
                if (inlig_pdb.split(".")[0] + "_min.pdb" not in os.listdir(datadir)):
                    get_Co(datadir,fn, inlig_pdb, st)
                else:
                    print("Use pervious generated Co")         
        print("Finish Co")
    
    return None
    
def feature_calculation_decoy(datadir, datadir_decoy, fn, ref_ligand_rdkit,ref_ligand_pdb,decoy_rdkit_list,decoy_pdb_list,water_type = "n"):
    '''
    feature calculation for decoys

    '''

    if water_type == "rw" or "w":
        ref_protein = fn + "_proten_RW.pdb"
        d_type = "_RW"
    else:
        ref_protein = fn + "_protein.pdb"
        d_type = ""
    ### copy ref protein ###
    cmd = "copy " + os.path.join(datadir,ref_protein) + " " + datadir_decoy
    os.system(cmd)
    

    outfile_V58= open(os.path.join(datadir_decoy,"Vina58_decoys" + d_type + ".csv"),"w")
    outfile_V58.write('pdb,idx,vina,' + ','.join(['vina' + str(n+1) for n in range(58)]) + "\n")
    outfile_SASA = open(os.path.join(datadir_decoy,"SASA_decoy" + d_type + ".csv"),"w")
    f_type = ["P","N","DA","D","A","AR","H","PL","HA","SA"]
    f_feature = ["P2." + i for i in f_type] + ["P2dl." + i for i in f_type] + ["P2dp." + i for i in f_type]
    outfile_SASA.write("pdb,idx," + ",".join(f_feature) + "\n")
    outfile_Ion = open(os.path.join(datadir_decoy,"Num_Ions" + d_type + ".csv"),"w")
    outfile_Ion.write("pdb,idx,Ni\n")

    for idx, decoy in enumerate(decoy_pdb_list):

        ### get Vina58 ###
        outfile = open(os.path.join(datadir_decoy,"Vina58" + str(idx+1) + ".csv"),"w")
        featureVina(outfile, fn, ref_protein, decoy, datadir_decoy)
        outfile.close()
        lines = open(os.path.join(datadir_decoy,"Vina58" + str(idx+1) + ".csv")).readlines()
        outfile_V58.write(fn + "," + str(idx + 1) + "," + ",".join(lines[0].split(",")[1:]))
        rm_cmd = "rm " + os.path.join(datadir_decoy,"Vina58" + str(idx+1) + ".csv")
        os.system(rm_cmd)
        print("Finish Vina" + str(idx+1))

        ### get SASA ###
        outfile  = open(os.path.join(datadir_decoy, "SASA" + str(idx+1) + ".csv"),"w")
        cal_SASA(outfile,fn,decoy,ref_protein,datadir_decoy)
        outfile.close()
        lines = open(os.path.join(datadir_decoy, "SASA" + str(idx+1) + ".csv")).readlines()
        outfile_SASA.write(fn + "," + str(idx + 1) + "," + ",".join(lines[0].split(",")[1:]))
        rm_cmd = "rm " + os.path.join(datadir_decoy,"SASA" + str(idx+1) + ".csv")
        os.system(rm_cmd)
        print("Finish SASA" + str(idx+1))

        ### get Ion ###
        outfile = open(os.path.join(datadir_decoy,"Num_Ions" + str(idx+1) + ".csv"),"w")
        cal_Ni(outfile,fn, ref_protein, decoy, datadir_decoy)
        outfile.close()
        lines = open(os.path.join(datadir_decoy,"Num_Ions" + str(idx+1) + ".csv")).readlines()
        outfile_Ion.write(fn + "," + str(idx + 1) + "," + ",".join(lines[0].split(",")[1:]))
        rm_cmd = "rm " + os.path.join(datadir_decoy,"Num_Ions" + str(idx+1) + ".csv")
        os.system(rm_cmd)
        print("Finish Ion" + str(idx+1))
    outfile_V58.close()
    outfile_SASA.close()
    outfile_Ion.close()


    ### get dERMSD ###
    outfile_dE = open(os.path.join(datadir_decoy,"dE_RMSD.csv"),"w")
    outfile_dE.write("pdb,idx,dE_global,RMSD_global,number0,number1\n")
    ### copy previous generated confs ###
    confs = os.path.join(datadir, fn + "_ligand_confs.sdf")
    lowest = os.path.join(datadir, fn + "_ligand_global_min.sdf")
    cmd = "cp " + confs + " " + datadir_decoy
    os.system(cmd)
    cmd = "cp " + lowest + " " + datadir_decoy
    for idx, decoy in enumerate(decoy_rdkit_list):
        outfile = open(os.path.join(datadir_decoy, "dE_RMSD" + str(idx + 1) + ".csv"),"w")
        feature_cal(outfile,fn, decoy, datadir_decoy, calc_type = "none", rewrite = False)
        outfile.close()
        lines = open(os.path.join(datadir_decoy, "dE_RMSD" + str(idx + 1) + ".csv")).readlines()
        outfile_dE.write(fn + "," + str(idx + 1) + "," + ",".join(lines[0].split(",")) + "\n")
        rm_cmd = "rm " + os.path.join(datadir_decoy,"dE_RMSD" + str(idx+1) + ".csv")
        os.system(rm_cmd)
        print("Finish dE_RMSD" + str(idx+1))
    outfile_dE.close()

    ### run fragments ###
    ### copy ref_lig ###
    cmd = "cp " + os.path(datadir, ref_ligand_pdb) +  " " + datadir_decoy
    os.system(cmd)
    cmd = "cp " + os.path(datadir, ref_ligand_rdkit) +  " " + datadir_decoy
    os.system(cmd)
    run_fragments(fn, datadir_decoy, ref_ligand_rdkit, ref_ligand_pdb, opt = False, water = False, decoy = True, decoy_list = decoy_pdb_list, decoy_type = d_type, decoy_pro = ref_protein)

    ### combine data ###
    combine(datadir,d_type,decoy = True)

    return None

def feature_calculation_ligand(datadir,fn, inlig_pdb, inlig_rdkit, inpro_pro, water_type, opt_type):
    '''
    feature calculation for ligands (C, Co, Crwo)
    
    '''

    ### update input structures ###
    inlig_C = inlig_pdb
    inlig_Co = fn + "_lig_min.pdb"
    inlig_Crwo = fn + "_lig_min_RW.pdb"

    inpro_C = inpro_pro
    inpro_Crw = fn + "_protein_RW.pdb"

    ################################

    if opt_type == "wo":
        d_type = {"":[inpro_C, inlig_C],"_min": [inpro_C, inlig_Co],"_min_RW": [inpro_Crw, inlig_Crwo]}
        ### get bw ###
        out_total = open(os.path.join(datadir,"Feature_BW_min_RW.csv"),"w")
        out_total.write("pdb,Nbw,Epw,Elw\n")
        cal_BW(out_total,fn,inpro_pro,inlig_Crwo,inpro_Crw,datadir)
        out_total.close()
        print("Finish BW")
    elif opt_type == "o":
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
    if opt_type == "wo":
        run_fragments(fn, datadir, inlig_rdkit, inlig_pdb,  opt = True, water = True)
    elif opt_type == "o"
        run_fragments(fn,datadir,inlig_rdkit, opt = True, water = False)
    else:
        run_fragments(fn,datadir,inlig_rdkit, opt = False, water = False)


    
    ### combine data ###
    for i in d_type.keys():
        combine(datadir,i)
    
    print("Finish Feature Calculation")

    return None







def run_features(datadir,fn, water_type = "rw", opt_type = "wo", ligand_type = "wo", rewrite = False, decoy = "CASF-docking"):
    '''

    run features
    
    water_type: receptor water type, defaults to "rw"
            "rw" --> get receptor water based our criteria
            "w"  --> all waters in protein_all.pdb are considered as receptor water
            "n"  --> no consideration of water molecules
    
    opt_type: optimization type, defaults to "wo"
            "wo" --> Crwo, Co
            "o"  --> Co
            "n"  --> no optimization 
            
    rewrite: whether to rewrite all features

    decoy: CASF decoys or not

    '''
    if decoy:
        ### CASF-2013/2016 docking/screening, no water has been used in decoy preparation ###
        opt_type = "n"
        water_type = "n"
        datadir_decoy = os.path.join(datadir,fn,"decoys")
        os.system("mkdir " + datadir_decoy)
        ref_ligand_rdkit, ref_ligand_pdb, decoy_rdkit_list, decoy_pdb_list = get_input_decoy(datadir,datadir_decoy, fm)
            
    else:
        inlig_rdkit, inlig_pdb, inpro_pro, inpro_water = get_input(datadir,fn)
        print("Finish Input Preparation")

    ### receptor water ###
    if water_type != "n":
        prepare_rw_receptor(datadir, fn, inpro_pro, inpro_water, water_type, rewrite)
    else:
        print("No Consideration of Water")

    ### get Co, Crwo ###
    if opt_type != "n":
        prepare_opt_ligand(datadir,fn,inlig_pdb,opt_type,rewrite)
    else:
        print("No Optimized Ligand")


    ### update input structures ###
    if decoy:
        feature_calculation_decoy(datadir, datadir_decoy, fn, ref_ligand_rdkit,ref_ligand_pdb,decoy_rdkit_list,decoy_pdb_list,water_type = water_type)
    else:
        feature_calculation_ligand(datadir,fn, inlig_pdb, inlig_rdkit, inpro_pro, water_type, opt_type)
    
    print("Finish Feature Calculation")

    return None


if __name__ == "__main__":
    datadir = "/Users/jianinglu1/Documents/GitHub/deltaVinaXGB_develop/Test"
    fn = "01"
    run_features(datadir,fn)
 