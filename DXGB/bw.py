#-----------------------------------------------------------------------------
# Bridging Water Features 
#-----------------------------------------------------------------------------
import os
import sys
import numpy as np
from DXGB.get_pdbinfo import *


#get the structural information of water molecules and get the bridging water molecules based on structure information
#output: bridging_water_info.txt BridgingWater_n.pdb

def get_angle(a,b,c):
    ba = a - b
    bc = c - b
    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.degrees(np.arccos(cosine_angle))

    return angle

def get_BW(fn,water,lig):
    outfile = open("BW_info.txt","w")

    lig = pdbinfo(fn,file = lig)
    lig = lig.getPolarAtoms()
    pro = pdbinfo(fn,file = water)
    pro_atoms = pro.getPolarAtoms()
    protein,waters = pdbinfo(fn, lines = pro_atoms).getProteinWaters()

    waters_coord_init = pdbinfo(fn,lines = waters).getCoords()
    protein_coord_init = pdbinfo(fn,lines = protein).getCoords()
    ligand_coord_init = pdbinfo(fn,lines=lig).getCoords()

    ### select the HOH in 3.5 of ligand/protein
    waters_coord = np.expand_dims(waters_coord_init, 1)
    protein_coord = np.expand_dims(protein_coord_init,0)
    ligand_coord = np.expand_dims(ligand_coord_init,0)

    if waters_coord.shape[0] == 0:
        print("No Bridging Water")
        outfile.close()
    else:
        distance_pw = np.linalg.norm(waters_coord - protein_coord, axis = 2)
        distance_lw = np.linalg.norm(waters_coord - ligand_coord, axis = 2)
        distance_pw_min = np.min(distance_pw, axis = 1)
        distance_lw_min = np.min(distance_lw, axis = 1)

        bw_index = []
        for idx, (i1,i2) in enumerate(zip(distance_pw_min,distance_lw_min)):
            if i1 > 2.0 and i1 < 3.5 and i2 > 2.0 and i2 < 3.5:
                bw_index.append(idx)
    
        for idx in bw_index:
            bw_coord = waters_coord_init[idx]
            bw_distance_p = distance_pw[idx]
            bw_distance_l = distance_lw[idx]
            for num1, i1 in enumerate(bw_distance_p):
                if i1 > 2.0 and i1 < 3.5:
                    p_coord = protein_coord_init[num1]
                    for num2, i2 in enumerate(bw_distance_l):
                        if i2 > 2.0 and i2 < 3.5:
                            l_coord = ligand_coord_init[num2]
                            angle = get_angle(p_coord, bw_coord, l_coord)
                            if angle >= 60:
                                bw_line = waters[idx]
                                bw_chain = chid(bw_line)
                                if bw_chain != " ":
                                    bw_name = str(int(resi(bw_line))) + "." + chid(bw_line)
                                else:
                                    bw_name = str(int(resi(bw_line)))
                                pro_line = protein[num1]
                                lig_line = lig[num2]
                                pro_chain = chid(pro_line)
                                pro_name = resn(pro_line)
                                if pro_chain != " ":
                                    pro_idx = str(int(resi(pro_line))) + "." + pro_chain 
                                else:
                                    pro_idx = str(int(resi(pro_line)))
                                pro_aname = atmn(pro_line).strip()
                                lig_name = atmn(lig_line).strip()
                                outline = fn + "," + pro_name + "," + pro_idx + "," + pro_aname + "," + bw_name + "," + lig_name + "," + str(round(i1,2)) + "," + str(round(i2,2)) + "," + str(int(round(angle))) + "\n"
                                outfile.write(outline)
        outfile.close()

    return None

def get_water(fn,water):
    #Get water residue index and water molecule file
    Residue_all = set([line.split(",")[4] for line in open("BW_info.txt")])
    print("BW satisfiles structural requirement:" + str(len(Residue_all)))
    index = open("water_index.txt","w")
    for i in set(Residue_all):
        index.write(i + "\n")
        inputwater = open(water)
        if "." in i:
            reid = i.split(".")[0]
            recid = i.split(".")[1]
            outfile = open("BW_" + reid + "_" + recid + ".pdb","w")
        else:
            reid = i
            recid = " "
            outfile = open("BW_" + reid + "_chain.pdb","w")
        for line_water in inputwater:
            if line_water[0:6] in ["HETATM","ATOM  "] and line_water[17:20] in ["WAT", "HOH"] and int(line_water[22:26])==int(reid) and line_water[21:22] == recid:
                outfile.write(line_water)
        outfile.write("END \n")
        outfile.close()
    index.close()

    return None

def addH(fn):
    ''' Add H to water molecule file (Vina need) '''
    for filename in os.listdir("."):
        if filename.startswith("BW") and filename.endswith(".pdb"):
            atoms = pdbinfo(file = filename).getAtoms()
            if len(atoms) == 3:
                continue
            else:
                newfilename = filename.split(".")[0] + "_addh.pdb"
                os.system("obabel " + filename + " -O " + newfilename + " -h")

                lines = [line for line in open(newfilename) if line[0:6] in ["ATOM  ", "HETATM"]]
                if "ATOM  " in lines[0]:
                    header = "ATOM  "
                else:
                    header = "HETATM"
                for i in range(1,3):
                    lines[i] = header + lines[i][6:]

                out = open(filename.split(".")[0] + "_addh_correct.pdb","w")
                out.write("".join(lines))
                out.close()

                os.system("mv " + filename.split(".")[0] + "_addh_correct.pdb" + " " + filename)
                os.system("rm " + newfilename)

    return None

#calculate vinascore
def genPDBQT(fn,pro,lig):
    olddir = os.getcwd()
    os.system("mkdir vina_BW")
    os.chdir("vina_BW")
    propdbqt = fn + "_prot_rc.pdbqt"
    ligpdbqt = fn + "_lig_rc.pdbqt"
    cmd1 ="$MGLPY $MGLUTIL/prepare_receptor4.py -r "  + pro + " -o " + propdbqt + " -U 'nphs' > out_lig.log"
    cmd2 ="$MGLPY $MGLUTIL/prepare_receptor4.py -r " + lig  + " -o " + ligpdbqt +  " -U 'nphs' > out_pro.log"
    os.system(cmd1)
    os.system(cmd2)
    for n, line in enumerate(open("../water_index.txt")):
        line = line.rstrip()
        if "." in line:
            wpdb = "../BW_" + line.split(".")[0] + "_" + line.split(".")[1] + ".pdb"
        else:
            wpdb = "../BW_" + line.split(".")[0] + "_chain.pdb"
        wpdbqt = "BW_" + str(n) + ".pdbqt"
        cmd3 ="$MGLPY $MGLUTIL/prepare_ligand4.py -l " + wpdb  + " -o " + wpdbqt +  " -U 'nphs' > out_bw_" + str(n) + ".log"
        cmd_pw = "$VINADIR/vina --receptor " + propdbqt + " --ligand " + wpdbqt + "  --score_only --log score_PW_" + str(n) + ".txt > out_pw_" + str(n) + ".log"
        cmd_lw = "$VINADIR/vina --receptor " + ligpdbqt + " --ligand " + wpdbqt + "  --score_only --log score_LW_" + str(n) + ".txt > out_lw_" + str(n) + ".log"
        try:
            os.system(cmd3)
            os.system(cmd_pw)
            os.system(cmd_lw)
        except:
            os.system("touch FAIL.log")
    os.chdir(olddir)

    return None

#get vina score information
#outfile: PWE_seperate.dat, LWE_seperate.dat, bridgingwater_vinascore_seperate.dat, bridgingwater_vinascore_total.dat
def get_result_PW(fn,out_PW):
    num = 0
    olddir = os.getcwd()
    os.chdir("vina_BW")
    num = len([ file for file in os.listdir(".") if file.startswith(("score_PW_")) and file.endswith(".txt")])
    for i in range(num):
        for line in open("score_PW_" + str(i) + ".txt"):
            if line[0:8] == "Affinity":
                try:
                    ### linux vina ###
                    value = float(line.split()[1])
                except:
                    ### os vina ###
                    value = float(line.split()[2])
                out_PW.write(fn + "," + str(i) + "," + str(value)  + "\n")
    os.chdir(olddir)
    out_PW.close()
    return None

def get_result_LW(fn,out_LW):
    num = 0
    olddir = os.getcwd()
    os.chdir("vina_BW")
    num = len([ file for file in os.listdir(".") if file.startswith(("score_LW_")) and file.endswith(".txt")])
    for i in range(num):
        for line in open("score_LW_" + str(i) + ".txt"):
            if line[0:8] == "Affinity":
                try:
                    ### linux vina ###
                    value = float(line.split()[1])
                except:
                    ### os vina ###
                    value = float(line.split()[2])
                out_LW.write(fn + "," + str(i) + "," + str(value)  + "\n")
    os.chdir(olddir)
    out_LW.close()
    return None

def get_BW_final(fn,out,out_total):
    value_PW = 0.00
    value_LW = 0.00
    index_list = []
    m = 0
    for line in open("Epw.csv"):
        if line.split(",")[0] == fn:
            if float(line.split(",")[2].strip("\n")) < 0:
                index = line.split(",")[1]
                for line_LW in open("Elw.csv"):
                    if line_LW.split(",")[0] == fn and line_LW.split(",")[1] == index and float(line_LW.split(",")[2].strip("\n")) < 0:
                        index_list.append(index)
                        out.write(fn + "," + index + "," + line.split(",")[2].strip("\n") + "," + line_LW.split(",")[2].strip("\n") + "\n")
                        value_PW += float(line.split(",")[2].strip("\n"))*-0.73349
                        value_LW += float(line_LW.split(",")[2].strip("\n"))*-0.73349
                        m = m + 1
    out.close()
    out_total.write(fn  + "," + str(m) + "," + str(value_PW) + "," + str(value_LW) + "\n")
    return index_list

#get the bridging water molecules
#outfile: Bridging_water_total.pdb, fn_protein_SF_bridgingwater.pdb
def get_waterfile(fn,pro, index):
    out = open("BW_total.pdb","w")
    for idx, line in enumerate(open("water_index.txt")):
        if str(idx) in index:
            if "." in line:
                file_index = "BW_" + "_".join(line.rstrip().split(".")) + ".pdb"
            else:
                file_index = "BW_" + line.rstrip() + "_chain.pdb"
            lines = [ l for l in open(file_index) if (l[0:6] in ["ATOM  ","HETATM"]) and (l[0:3] != "END")]
            out.write("".join(lines))
    out.close()

    out = open("BW_total.pdb")
    count_lines = int(len(out.readlines())/3)
    out.close()

    if os.stat("BW_total.pdb").st_size == 0:
        print('No BW has been found in ' + fn)
    else:
        print(str(count_lines) +" BW have been saved in " + fn  + "_protein_BW.pdb")
    
    in_rec = open(pro)
    out_rec = open(fn + "_protein_BW.pdb","w")
    lines = [line for line in in_rec if line[0:6] in ["ATOM  ", "HETATM"]]
    out_rec.write("".join(lines))
    in_water = open("BW_total.pdb")
    lines = [line for line in in_water]
    out_rec.write("".join(lines))
    out_rec.write("END")
    out_rec.close()
    return None

def cal_BW(out_total,fn,inprot,inlig,inwater,datadir, Feature = True):
    os.chdir(datadir)
    inidir = os.getcwd()
    pro = os.path.join(datadir,inprot)
    water =os.path.join(datadir,inwater)
    lig = os.path.join(datadir,inlig)
    os.system("mkdir BW")
    os.chdir("BW")
    get_BW(fn,water,lig)
    get_water(fn,water)
    addH(fn)
    genPDBQT(fn,pro,lig)
    out_PW = open("Epw.csv","w")
    out_LW = open("Elw.csv","w")
    get_result_PW(fn,out_PW)
    get_result_LW(fn,out_LW)
    out = open("BW.csv","w")
    index = get_BW_final(fn,out,out_total)
    get_waterfile(fn,pro,index)
    if not Feature:
        os.system("cp " + fn + "_protein_BW.pdb ../")
    os.chdir(inidir)
    os.system("rm -r BW*")

if __name__ == "__main__":
    ### test ###
    fn = "3c2f"
    out_total = open("/Users/jianinglu1/Documents/script/deltaXGB_linux/Feature/test/Feature_BW.csv","w")
    datadir = "/Users/jianinglu1/Documents/script/deltaXGB_linux/Feature/test"
    inprot =  fn + "_protein.pdb"
    inlig =  fn + "_ligand.pdb"
    inwater = fn + "_protein_all.pdb"
    cal_BW(out_total,fn,inprot,inlig,inwater,datadir)
