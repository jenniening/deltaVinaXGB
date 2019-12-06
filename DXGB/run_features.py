import os
import sys

import rdkit 
from rdkit import Chem
import DXGB.get_pdbinfo as get_pdbinfo

import DXGB.bw
from DXGB.bw import cal_BW
import DXGB.rw 
from DXGB.rw import get_Crw
import DXGB.opt
from DXGB.opt import get_Co
import DXGB.cal_vina58
from DXGB.cal_vina58 import featureVina
import DXGB.cal_sasa
from DXGB.cal_sasa import cal_SASA
import DXGB.cal_dERMSD
from DXGB.cal_dERMSD import feature_cal
import DXGB.cal_ion
from DXGB.cal_ion import cal_Ni

import DXGB.combine_data
from DXGB.combine_data import combine


def renumber(fmt, infile, outfile):
    """
    Rename atoms in file based on order 
    """
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


def get_input(datadir, fn):
    """
    Get input files based on pdbid
    return: inlig_rdkit --> inlig mol2 or sdf for RDkit using
            inlig3 --> inlig pdb file
            inpro1 --> inpro with only protein
            inpro2 --> inpro with both proteina and water
    """
    olddir = os.getcwd()
    os.chdir(datadir)
    ### get ligand ###
    ### input should be provided as either of sdf file (best choice) or mol2 file ###
    inlig1 = fn + "_ligand.mol2"
    inlig2 = fn + "_ligand.sdf"
    inlig3 = fn + "_ligand.pdb"
    ### check ligand input file ###
    inlig_rdkit = None
    if inlig1 in os.listdir("."):
        inlig = inlig1
        try:
            mol = Chem.MolFromMol2File(inlig,removeHs = False)
            if mol == None:
                if inlig2 in os.listdir("."):
                    inlig = inlig2
                    mol = Chem.SDMolSupplier(inlig,removeHs = False)[0]
                    if mol != None:
                        inlig_rdkit = inlig2
            else:
                inlig_rdkit = inlig1
        except:
            pass
    else:
        inlig = inlig2
        try:
            mol = Chem.SDMolSupplier(inlig,removeHs = False)[0]
            if mol != None:
                inlig_rdkit = inlig2
        except:
            pass
    ### correcting atom name for sasa calculation ###
    ### if inlig sdf or mol2 file might have problems, you should also provide a pdb file that can be used to conduct other calculation ###
    if inlig3 not in os.listdir(".") and inlig_rdkit == None:
        ### if sdf and mol2 files can't be processed by RDKit, we need to generate pdb file and continue other calculations ###
        ### this is not good, since if that molecule has large problem, the conversion process might be wrong; but this can work when the problem is caused by RDKit ###
        infmt = inlig.split(".")[-1]
        outlig = inlig3
        cmd = "obabel -i" + infmt + " " +  outlig + " -opdb -O " + outlig
        os.system(cmd)

    if inlig3 not in os.listdir("."):
        inlig = inlig_rdkit
        infmt = inlig.split(".")[-1]
        if infmt == "mol2":
            outlig_num = inlig.split(".")[0] + "_rename.mol2"
            renumber(infmt,inlig,outlig_num)
            outlig = inlig3.split(".")[0] + "_rename.pdb"
            cmd = "obabel -i" + infmt + " " +  outlig_num + " -opdb -O " + outlig
            os.system(cmd)
        elif infmt == "sdf":
            outlig = inlig3
            cmd = "obabel -i" + infmt + " " +  inlig + " -opdb -O " + outlig
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
    if inlig_rdkit != None and os.path.isfile(inlig_rdkit) and os.stat(inlig_rdkit).st_size != 0:
        print("Ligand for conformation stability:" + inlig_rdkit)
    else:
        print("Warning:input ligand should be checked, skip ligand stability calculation, use default(dE:-300, RMSD:300)")

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


def prepare_rw_receptor(datadir, fn, inpro_pro, inpro_water, inlig, water_type, rewrite = False):
    """
    prepare protein with water for different water type
    :param datadir: directory for files
    :param fn: index for input 
    :param inpro_pro: inpro only protein file 
    :param inpro_water: inpro with protein and water molecules file
    :param inlig: inlig pdb file 
    :param water_type: protein part water type, defaults to "rbw"
            "rbw" --> get protein part water based RW and BW
            "rw" --> get protein part water based RW
            "bw" --> get protein part water based BW
            "pw" --> get protein part water based PW (all waters in protein_all.pdb are considered as protein part water)
            "n"  --> no consideration of water molecules
    :param rewrite: whether to rewrite RW structure generated previously, defaults to False

    """
    
    if water_type == "rbw":
        print("Protein Water: calculate both RW and BW")
        if rewrite:
            get_Crw(fn,inpro_pro,inpro_water,datadir)
            print("Finish generate RW")
            out_total = open("Feature_BW_initial.csv","w")
            cal_BW(out_total,fn,inpro_pro,inlig,inpro_water,datadir, Feature = False)
            out_total.close()
            print("Finish generate BW ")
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
        print("Protein Water: calculate by RW")
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
        print("Protein Water: calculate BW")
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


def prepare_opt(datadir, fn, inlig_pdb, opt_type, rewrite = False):
    """
    Prepare AutoDock Vina optimized ligand for different optimization type

    :param datadir: directory for files
    :param fn: index for input 
    :param inlig_pdb: inlig pdb file
    :param opt_type: based on receptor water type, defaults to "rbwo"
            "rbwo" --> do optimization for protein without water, with RW and with BW
            "rwo" --> do optimization for protein with RW
            "bwo" --> do optimization for protein with BW
            "pwo" --> do optimization for protein with PW(all waters in protein_all.pdb are considered as protein part water)
            "o" --> do optimization for protein without water
            "n" --> no optimization
    :param rewrite: whether to rewrite optimized structure generated previously, defaults to False

    """
    if opt_type == "rbwo":
        for st in ["","_RW","_BW"]:
            if rewrite:
                get_Co(datadir,fn, inlig_pdb, st)
            else:
                if fn + "_lig_min" + st + ".pdb" not in os.listdir(datadir):
                    get_Co(datadir,fn, inlig_pdb, st)
                else:
                    print("Use pervious generated C" + st + "O")
        print("Finish Optimization")

    elif opt_type == "rwo":
        for st in ["_RW"]:
            if rewrite:
                get_Co(datadir,fn, inlig_pdb, st)
            else:
                if fn + "_lig_min" + st + ".pdb" not in os.listdir(datadir):
                    get_Co(datadir,fn, inlig_pdb, st)
                else:
                    print("Use pervious generated C" + st + "O")
        print("Finish Optimization")

    elif opt_type == "bwo":
        for st in ["_BW"]:
            if rewrite:
                get_Co(datadir,fn, inlig_pdb, st )
            else:
                if fn + "_lig_min" + st + ".pdb" not in os.listdir(datadir):
                    get_Co(datadir,fn, inlig_pdb, st )
                else:
                    print("Use pervious generated C" + st + "O")
        print("Finish Optimization")

    elif opt_type == "pwo":
        for st in ["_PW"]:
            if rewrite:
                get_Co(datadir,fn, inlig_pdb, st)
            else:
                if fn + "_lig_min" + st + ".pdb" not in os.listdir(datadir):
                    get_Co(datadir,fn, inlig_pdb, st)
                else:
                    print("Use pervious generated C" + st + "O")
        print("Finish Optimization")

    elif opt_type == "o":
        for st in [""]:
            if rewrite:
                get_Co(datadir,fn, inlig_pdb, st)
            else:
                if (inlig_pdb.split(".")[0] + "_min.pdb" not in os.listdir(datadir)):
                    get_Co(datadir,fn, inlig_pdb, st)
                else:
                    print("Use pervious generated C" + st + "O")         
        print("Finish Optimization")


def run_Vina_features(datadir, d_type, fn, inpro, inlig):
    """
    Calculate Vina features

    :param datadir: directory for input file
    :param d_type: structure type after consideration opt and water
    :param fn: input index
    :param inpro: protein file
    :param inlig: ligand file

    """
    outfile = open(os.path.join(datadir,"Vina58" + d_type + ".csv"),"w")
    outfile.write('pdb,vina,' + ','.join(['vina' + str(n+1) for n in range(58)]) + "\n")
    featureVina(outfile, fn, inpro, inlig, datadir)
    outfile.close()


def run_SASA_features(datadir, d_type, fn, inpro, inlig):
    """
    Calculate SASA features
    
    :param datadir: directory for input file
    :param d_type: structure type after consideration opt and water
    :param fn: input index
    :param inpro: protein file
    :param inlig: ligand file

    """
    f_type = ["P","N","DA","D","A","AR","H","PL","HA","SA"]
    f_feature = ["P2." + i for i in f_type] + ["P2dl." + i for i in f_type] + ["P2dp." + i for i in f_type]
    out_SASA = open(os.path.join(datadir, "SASA" + d_type + ".csv"),"w")
    out_SASA.write("pdb," + ",".join(f_feature) + "\n")
    cal_SASA(out_SASA,fn,inlig,inpro,datadir)
    out_SASA.close()


def run_Ion_features(datadir, d_type, fn, inpro, inlig):
    """
    Calculate Ion features

    :param datadir: directory for input file
    :param d_type: structure type after consideration opt and water
    :param fn: input index
    :param inpro: protein file
    :param inlig: ligand file

    """
    outfile = open(os.path.join(datadir,"Num_Ions" + d_type + ".csv"),"w")
    outfile.write("pdb,Ni\n")
    cal_Ni(outfile,fn, inpro, inlig, datadir)
    outfile.close()


def run_BW_features(datadir, d_type, fn, inpro_pro, inlig, inpro):
    """
    Calculate BW features 
    
    :param datadir: directory for input file
    :param d_type: structure type after consideration opt and water
    :param fn: input index
    :param inpro_pro: protein structure with only protein
    :param inlig: ligand file
    :param inpro: protein file

    """

    out_total = open(os.path.join(datadir,"Feature_BW" + d_type + ".csv"),"w")
    out_total.write("pdb,Nbw,Epw,Elw\n")
    cal_BW(out_total,fn,inpro_pro,inlig,inpro,datadir)
    out_total.close()


def run_dE_features(datadir, fn, inlig_rdkit, rewrite = False):
    """
    Calculate dE features

    :param datadir: directory for input file
    :param fn: input index
    :param inlig_rdkit: inlig mol2 or sdf file for RDkit
    :param rewrite: whether to rewrite conformation file generated previously, defaults to False

    """
    outfile = open(os.path.join(datadir,"dE_RMSD.csv"),"w")
    outfile.write("pdb,dE_global,RMSD_global\n")
    if inlig_rdkit:
        feature_cal(outfile,fn, inlig_rdkit, datadir, calc_type = "GenConfs", rewrite = rewrite)
    else:
        outfile.write(fn + ",-300,300\n")
    outfile.close()


def feature_calculation_ligand(datadir,fn, inlig_pdb, inlig_rdkit, inpro_pro, water_type, opt_type, rewrite = False, feature_type = "all"):
    """
    Feature calculation for ligands (C, Co, Crwo, Cbwo, Cpwo)
    
    :param datadir: datadir for structure
    :param fn: input index
    :param inlig_pdb: inlig file ends with .pdb
    :param inlig_rdkit: inlig file for ligand stability that is calculated using RDkit
    :param inpro_pro: inpro file 
    :param water_type: water type 
    :param opt_type: opt type
    :param rewrite: whether to rewrite generated conformations
    :param feature_type: which feature will be calculated, defaults to "all"
    :return: [description]
    :rtype: [type]

    """
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
        if water_type == "rbw":
            d_type = {"":[inpro_C, inlig_C], "_RW":[inpro_Crw, inlig_C], "_BW":[inpro_Cbw, inlig_C]}
        elif water_type == "rw":
            d_type = {"_RW":[inpro_Crw, inlig_C]}
        elif water_type == "bw":
            d_type = {"_BW":[inpro_Cbw, inlig_C]}
        elif water_type == "pw":
            d_type = {"_PW":[inpro_Cpw, inlig_C]}
        else:
            d_type = {"":[inpro_C, inlig_C]}

    for i in d_type.keys():
        if i == "":
            print("C")
        elif i == "_min":
            print("Co")
        elif i == "_RW":
            print("Crw")
        elif i == "_min_RW":
            print("Crwo")
        elif i == "_BW":
            print("Cbw")
        elif i == "_min_BW":
            print("Cbwo")
        elif i == "_PW":
            print("Cpw")
        elif i == "_min_PW":
            print("Cpwo")
        inpro = d_type[i][0]
        inlig = d_type[i][1]
        ### get BW features ###
        ### only for structure with water in protein ###
        if i in ["_min_RW", "_min_BW", "_min_PW", "_RW","_BW","_PW"]:
            inpro_pro = inpro_C
            if feature_type == "all" or feature_type == "BW":
                run_BW_features(datadir, i, fn, inpro_pro, inlig, inpro)
                print("Finish Bridging Water feature calculation, save in Feature_BW" + i + ".csv")
            else:
                print("Use previous calculated Bridging Water feature in Feature_BW" + i + ".csv")
        ### get Vina58 ###
        if feature_type == "all" or feature_type == "Vina":
            run_Vina_features(datadir, i, fn, inpro, inlig)
            print("Finish Vina, save in Vina58" + i + ".csv")
        else:
            print("Use previous calculated Vina in Vina58" + i + ".csv")
        ### get sasa ###
        if feature_type == "all" or feature_type == "SASA":
            run_SASA_features(datadir, i, fn, inpro, inlig)
            print("Finish SASA, save in SASA" + i + ".csv")
        else:
            print("Use previous calculated SASA in SASA" + i + ".csv")
        ### get ion ###
        if feature_type == "all" or feature_type == "Ion":
            run_Ion_features(datadir, i, fn, inpro, inlig)
            print("Finish Ion, save in Num_Ions" + i + ".csv")
        else:
            print("Use previous calculated Ion in Num_Ions" + i + ".csv")

    ### get dERMSD ###
    if feature_type == "all" or feature_type == "dE":
        run_dE_features(datadir, fn, inlig_rdkit, rewrite)
        print("Finish ligand stability calculation, save in dE_RMSD.csv")
    else:
        print("Use previous calculated ligand stability in dE_RMSD.csv")

    ### combine data ###
    if feature_type == "all":
        for i in d_type.keys():
            combine(datadir,i)
            print("Features has been saved in Input" + i + ".csv")


def run_features(datadir, pdbid, water_type = "rbw", opt_type = "rbwo", rewrite = False, feature_type = "all"):
    """
    Calculate features
    
    datadir: directory for structures
    pdbid: input index

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
  
    rewrite: whether to rewrite structures generated for RW or conformations

    feature_type: which feature will be calcuated, defaults to "all"
            "Vina" --> only calculate Vina58 features
            "SASA" --> only calculate SASA features
            "BW" --> only calculate BW features
            "dE" --> only calculate ligand satbility features
            "Ion" --> only calculate Ion features
    
    """

    inlig_rdkit, inlig_pdb, inpro_pro, inpro_water = get_input(datadir,pdbid)
    print("Finish Input Preparation")

    if water_type != "n":
        prepare_rw_receptor(datadir, pdbid, inpro_pro, inpro_water, inlig_pdb, water_type, rewrite = rewrite)
        print("Consideration of Water Effect")
    else:
        print("No Consideration of Water")

    ### get Co, Crwo ###
    if opt_type != "n":
        prepare_opt(datadir, pdbid, inlig_pdb, opt_type, rewrite = rewrite)
    else:
        print("No Optimized Ligand")

    ### update input structures ###
    feature_calculation_ligand(datadir, pdbid, inlig_pdb, inlig_rdkit, inpro_pro, water_type, opt_type, rewrite, feature_type)
    print("Finish Feature Calculation")

    
 
