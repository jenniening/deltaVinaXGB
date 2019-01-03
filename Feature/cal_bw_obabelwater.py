"""
Get bridging water features
"""
__author__ = "Jianing Lu"
__copyright__ = "Copyright 2018, NYU"
__license__ = ""

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

from chimera import selection
import os
import chimera
from chimera import runCommand as rc
import sys
from software_path import path_mgl_python
from software_path import path_mgl_script
from software_path import path_vina
from software_path import path_obabel

MGLPY = path_mgl_python()

MGLUTIL = path_mgl_script()

vina = path_vina()
obable = path_obabel()

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------

""" Part1 get the structural information of water molecules and get the bridging water molecules based on structure information """
""" Outfile: bridging_water_info.txt BridgingWater_n.pdb """

def get_BW(fn,water,lig,pro):
    outfile = open("BW_info.txt","w")
    w = chimera.openModels.open(water)[0]
    l = chimera.openModels.open(lig)[0]
    p = chimera.openModels.open(pro)[0]
    # select the HOH in 3.5 of ligand/protein
    for line in open(water):
        if "HOH" in line:
            water_name = "HOH"
            break
        elif "WAT" in line:
            water_name = "WAT"
            break
    rc("sel #1@N=,O=,S= za<3.5 & #1@N=,O=,S= za>2.0 & #2@N=,O=,S= za>2.0 & #2@N=,O=,S= za<3.5 & #0:" + water_name + "@O=")
    #rc("sel #1@N=,O=,S= za<3.5 & #1@N=,O=,S= za>2.0 & #2@N=,O=,S= za>2.0 & #2@N=,O=,S= za<3.5 & #0:HOH@O=")
    wa =  selection.currentAtoms()
    # select the ligand in 3.5 of HOH
    rc("sel sel za<3.5 & #1@N=,O=,S=")
    la = selection.currentAtoms()
    # select the protein in 3.5 of HOH
    ### for oabebl convert water
    rc("sel #1@N=,O=,S= za<3.5 & #1@N=,O=,S= za>2.0 & #2@N=,O=,S= za>2.0 & #2@N=,O=,S= za<3.5 & #0:" + water_name + "@O=")
    #rc("sel #1@N=,O=,S= za<3.5 & #1@N=,O=,S= za>2.0 & #2@N=,O=,S= za>2.0 & #2@N=,O=,S= za<3.5 & #0:HOH@O=")
    rc("sel sel za<3.5 & #2@N=,O=,S=")
    pa = selection.currentAtoms()

    for a1 in pa:
        for a2 in wa:
            for a3 in la:
                d1 = a1.xformCoord().distance(a2.xformCoord())
                d2 = a2.xformCoord().distance(a3.xformCoord())
                if d1 <= 3.5 and d1 >= 2.0 and d2 <= 3.5 and d2 >= 2.0:
                    angle = chimera.angle(a1.xformCoord(), a2.xformCoord(), a3.xformCoord())
                    if angle >= 60:
                        rc("sel " + str(a1.residue))
                        tmp = selection.currentResidues()
                        for r in tmp:
                            prn = r.type
                            out = [prn, str(a1.residue)[3:],  a1.name, str(a2.residue)[3:], a3.name, str(round(d1,2)), str(round(d2,2)), str(int(angle))]
                            out =  ','.join(out)
                            outfile.write(fn + "," + out + "\n")
    outfile.close()
    rc("del")

def get_water(fn,water):
    inputfile = open("BW_info.txt")
    Residue_all = []
    for line in inputfile:
        residue = (line.split(",")[4]).split(".")[0]
        Residue_all.append(residue)
    index = open("water_index.txt","w")
    for i in sorted(set(Residue_all),key = int):
        index.write(i + "\n")
        inputwater = open(water)
        outfile = open("BW_" + i + ".pdb","w")
        for line_water in inputwater:
            ### for obabel convert water
            if line_water[0:6] in ["HETATM","ATOM  "]  and line_water[17:20] in ["WAT","HOH"] and int(line_water[22:26])==int(i):
            #if line_water[0:6] == "HETATM" and line_water[17:20] == "HOH" and int(line_water[22:26])==int(i):
                outfile.write(line_water)
        outfile.write("END \n")
        outfile.close()
    index.close()
    return None


def addH(fn):
    for filename in os.listdir("."):
        if filename.startswith("BW") and filename.endswith(".pdb"):
            rc("open " + filename)
            rc("addh")
            rc("write format pdb #0 " + filename)
            rc("del")
    return None


#####Part3 calculate the vinascore
def genPDBQT(fn,pro,lig):
    olddir = os.getcwd()
    os.system("mkdir vina_BW")
    os.chdir("vina_BW")
    propdbqt = fn + "_prot_rc.pdbqt"
    ligpdbqt = fn + "_lig_rc.pdbqt"
    cmd1 =MGLPY + " " + MGLUTIL +  "prepare_receptor4.py -r "  + pro + " -o " + propdbqt + " -U 'nphs' > out_lig.log"
    cmd2 =MGLPY + " " + MGLUTIL + "prepare_receptor4.py -r " + lig  + " -o " + ligpdbqt +  " -U 'nphs' > out_pro.log"
    os.system(cmd1)
    os.system(cmd2)
    n = 0
    for line in open("../water_index.txt"):
        wpdb = "../BW_" + line.strip("\n") + ".pdb"
        wpdbqt = "BW_" + str(n) + ".pdbqt"
        cmd3 =MGLPY + " " + MGLUTIL + "prepare_ligand4.py -l " + wpdb  + " -o " + wpdbqt +  " -U 'nphs' > out_bw_" + str(n) + ".log"
        cmd_pw =vina + " --receptor " + propdbqt + " --ligand " + wpdbqt + "  --score_only --log score_PW_" + str(n) + ".txt > out_pw_" + str(n) + ".log"
        cmd_lw =vina + " --receptor " + ligpdbqt + " --ligand " + wpdbqt + "  --score_only --log score_LW_" + str(n) + ".txt > out_lw_" + str(n) + ".log"
        try:
            os.system(cmd3)
            os.system(cmd_pw)
            os.system(cmd_lw)
        except:
            os.system("touch FAIL.log")
        n = n + 1
    os.chdir(olddir)
    return None



#####Part4 Get the vina score information
#####outfile: PWE_seperate.dat, LWE_seperate.dat, bridgingwater_vinascore_seperate.dat, bridgingwater_vinascore_total.dat
def get_result_PW(fn,out_PW):
    num = 0
    olddir = os.getcwd()
    os.chdir("vina_BW")
    for file in os.listdir("."):
        if file.startswith(("score_PW_")) and file.endswith(".txt"):
            num = num + 1
    for i in range(num):
        for line in open("score_PW_" + str(i) + ".txt"):
            if line[0:8] == "Affinity":
                ### linux vina
                value = float(line.split()[1])
                ###
                #value = float(line.split()[2])
                out_PW.write(fn + "," + str(i) + "," + str(value)  + "\n")
    os.chdir(olddir)
    out_PW.close()
    return None


def get_result_LW(fn,out_LW):
    num = 0
    olddir = os.getcwd()
    os.chdir("vina_BW")
    for file in os.listdir("."):
        if file.startswith(("score_LW_")) and file.endswith(".txt"):
            num = num +1
    for i in range(num):
        for line in open("score_LW_" + str(i) + ".txt"):
            if line[0:8] == "Affinity":
                ### linux vina
                value = float(line.split()[1])
                ###
                #value = float(line.split()[2])
                out_LW.write(fn + "," + str(i) + "," + str(value)  + "\n")
    os.chdir(olddir)
    out_LW.close()
    return None



def get_BW_final(fn,out,out_total):
    value_PW = 0.00
    value_LW = 0.00
    m = 0
    for line in open("Epw.csv"):
        if line.split(",")[0] == fn:
            if float(line.split(",")[2].strip("\n")) < 0:
                index = line.split(",")[1]
                for line_LW in open("Elw.csv"):
                    if line_LW.split(",")[0] == fn and line_LW.split(",")[1] == index and float(line_LW.split(",")[2].strip("\n")) < 0:
                        out.write(fn + "," + index + "," + line.split(",")[2].strip("\n") + "," + line_LW.split(",")[2].strip("\n") + "\n")
                        value_PW = value_PW + float(line.split(",")[2].strip("\n"))*-0.73349
                        value_LW = value_LW + float(line_LW.split(",")[2].strip("\n"))*-0.73349
                        m = m + 1
    out.close()
    out_total.write(fn  + "," + str(m) + "," + str(value_PW) + "," + str(value_LW) + "\n")
    return None

#####Part5 based on vina score information to get the bridging water molecules
#####outfile: Bridging_water_total.pdb, fn_protein_SF_bridgingwater.pdb
def get_waterindex():
    index = []
    for line in open("BW.csv"):
        index.append(line.split(",")[1])
    return index


def get_fileindex():
    file_index = []
    for line in open("water_index.txt"):
        file = "BW_" + line.strip("\n") + ".pdb"
        file_index.append(file)
    return file_index


def get_waterfile(fn,pro, index, file_index):
    out = open("BW_total.pdb","w")
    for i in index:
        for line in open(file_index[int(i)]):
            if line[0:6] == "ATOM  " or line[0:6] == "HETATM":
                if line[0:3] != "END":
                    out.write(line)
    out.close()
    in_rec = open(pro)
    out_rec = open(fn + "_protein_BW.pdb","w")
    for line in in_rec:
        if line[0:6] == "ATOM  " or line[0:6] == "HETATM":
            out_rec.write(line)
    in_water = open("BW_total.pdb")
    for line in in_water:
        out_rec.write(line)
    out_rec.write("END")
    out_rec.close()
    return None



def cal_BW(out_total,fn,inprot,inlig,inwater,datadir):
    os.chdir(datadir)
    inidir = os.getcwd()
    pro = datadir + inprot
    water = datadir + inwater
    lig = datadir + inlig
    os.system("mkdir BW")
    os.chdir("BW")
    get_BW(fn,water,lig,pro)
    get_water(fn,water)
    addH(fn)
    genPDBQT(fn,pro,lig)
    out_PW = open("Epw.csv","w")
    out_LW = open("Elw.csv","w")
    get_result_PW(fn,out_PW)
    get_result_LW(fn,out_LW)
    out = open("BW.csv","w")
    get_BW_final(fn,out,out_total)
    index = get_waterindex()
    file_index = get_fileindex()
    get_waterfile(fn,pro,index,file_index)
    os.chdir(inidir)
    os.system("rm -r BW*")




def main():
    args = sys.argv[1:]
    if args[-1] == "file":
        pdbfile = open('%s'%(args[0] + args[1]),'r')
        pdblist = []
        for i in pdbfile.readlines():
            pdblist.append(i[0:4])
    else:
        pdblist = []
        pdblist.append(args[1])
    datadir = args[0]
    out_total_min = open(datadir + "Feature_BW_min_RW.csv","w")
    out_total_min.write("pdb,Nbw,Epw,Elw\n")
    for fn in pdblist:
        inpro = fn + "_protein.pdb"
        inpro_water = fn + "_protein_RW.pdb"
        if fn + "_ligand.mol2" in os.listdir(datadir):
            inlig = fn + "_ligand.mol2"
        elif fn + "_ligand.sdf" in os.listdir(datadir):
            inlig =  fn + "_ligand.sdf"
            inlig_out =  fn + "_ligand.mol2"
            os.system(obable + " -isdf " + datadir + inlig + " -omol2 -O " + datadir + inlig_out)
            inlig = inlig_out
        elif fn + "_ligand_rigid.pdbqt" in os.listdir(datadir):
            inlig =  fn + "_ligand_rigid.pdbqt"
            inlig_out =  fn + "_ligand.mol2"
            os.system(obable + " -ipdbqt " + datadir + inlig + " -omol2 -O " + datadir + inlig_out)
            inlig = inlig_out
        elif fn + "_ligand_flexible.pdbqt" in os.listdir(datadir):
            inlig =  fn + "_ligand_flexible.pdbqt"
            inlig_out =  fn + "_ligand.mol2"
            os.system(obable + " -ipdbqt " + datadir + inlig + " -omol2 -O " + datadir +  inlig_out)
            inlig = inlig_out

        else:
            print("wrong ligand input file format, it should be mol2, sdf, or pdbqt")
        inlig_min = fn + "_lig_min_RW.pdb"
        cal_BW(out_total_min,fn,inpro,inlig_min,inpro_water,datadir)
    out_total_min.close()

    return None

if __name__ == "__main__":
    main()
