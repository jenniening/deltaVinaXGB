#-----------------------------------------------------------------------------
# Crw
#-----------------------------------------------------------------------------
import numpy as np
import os
import sys
import DXGB.get_pdbinfo as get_pdbinfo

def get_RW(fn, inpro):
    ''' Select the HOH in [2.0, 3.5] of protein '''

    outfile = open("RW_info.txt","w")
    pro = get_pdbinfo.pdbinfo(fn,file = inpro)
    pro_atoms = pro.getPolarAtoms()
    protein,waters = get_pdbinfo.pdbinfo(fn, lines = pro_atoms).getProteinWaters()
    waters_coord = get_pdbinfo.pdbinfo(fn,lines = waters).getCoords()
    protein_coord = get_pdbinfo.pdbinfo(fn,lines = protein).getCoords()
    ### calculate distance ###
    waters_coord = np.expand_dims(waters_coord, 1)
    protein_coord = np.expand_dims(protein_coord,0)
    if waters_coord.shape[0] == 0:
        print("No Receptor Water")
        outfile.close()
    else:
        distance = np.linalg.norm(waters_coord - protein_coord, axis = 2)
        distance_min = np.min(distance, axis = 1)
        rw_index = []
        for idx, i in enumerate(distance_min):
            if i > 2.0 and i < 3.5:
                rw_index.append(idx)
        for i in rw_index:
            rw_line = waters[i]
            rw_distance = distance[i]
            rw_chain = get_pdbinfo.chid(rw_line)
            if rw_chain != " ":
                rw_name = str(int(get_pdbinfo.resi(rw_line))) + "." + get_pdbinfo.chid(rw_line)
            else:
                rw_name = str(int(get_pdbinfo.resi(rw_line)))
            
            for idx, d in enumerate(rw_distance):
                if d < 3.5 and d > 2.0:
                    pro_line = protein[idx]
                    pro_chain = get_pdbinfo.chid(pro_line)
                    pro_name = get_pdbinfo.resn(pro_line)
                    if pro_chain != " ":
                        pro_idx = str(int(get_pdbinfo.resi(pro_line))) + "." + pro_chain 
                    else:
                        pro_idx = str(int(get_pdbinfo.resi(pro_line)))
                    pro_aname = get_pdbinfo.atmn(pro_line).strip()
                    outline = fn + "," + pro_name + "," + pro_idx + "," + pro_aname + "," + rw_name + "," + str(round(d,2)) + "\n"
                    outfile.write(outline)
        outfile.close()

def get_water(fn,water):
    ''' Get water residue index and water molecule file '''
 
    Residue_all = set([line.split(",")[4] for line in open("RW_info.txt")])
    print("RW satisfiles distance requirement:" + str(len(Residue_all)))
    index = open("water_index.txt","w")
    for i in set(Residue_all):
        index.write(i + "\n")
        inputwater = open(water)
        if "." in i:
            reid = i.split(".")[0]
            recid = i.split(".")[1]
            outfile = open("RW_" + reid + "_" + recid + ".pdb","w")
        else:
            reid = i
   
            recid = " "
            outfile = open("RW_" + reid + "_chain.pdb","w")
        for line_water in inputwater:
            if line_water[0:6] in ["HETATM","ATOM  "] and line_water[17:20] in ["WAT", "HOH"] and int(line_water[22:26])==int(reid) and line_water[21:22] == recid:
                outfile.write(line_water)
        outfile.write("END \n")
        outfile.close()
    index.close()


def addH(fn):
    ''' Add H to water molecule file (Vina need) '''

    for filename in os.listdir("."):
        if filename.startswith("RW") and filename.endswith(".pdb"):
            atoms = get_pdbinfo.pdbinfo(file = filename).getAtoms()
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


def runVina(fn,pro):
    ''' Get pdbqt file of protein, water and run Vina '''

    olddir = os.getcwd()
    os.system("mkdir vina_RW")
    os.chdir("vina_RW")
    
    ### get pdbqt file for protein ###
    propdbqt = fn + "_prot_rc.pdbqt"
    cmd1 = "$MGLPY $MGLUTIL/prepare_receptor4.py -r "  + pro + " -o " + propdbqt + " -U 'nphs' > out_pro.log"
    os.system(cmd1)
    
    ### get pdbqt files for water molecules and run Vina ###
    for n, line in enumerate(open("../water_index.txt")):
        line = line.rstrip()
        if "." in line:
            wpdb = "../RW_" + line.split(".")[0] + "_" + line.split(".")[1] + ".pdb"
        else:
            wpdb = "../RW_" + line.split(".")[0] + "_chain.pdb"
        wpdbqt = "RW_" + str(n) + ".pdbqt"
        cmd2 = "$MGLPY $MGLUTIL/prepare_ligand4.py -l " + wpdb  + " -o " + wpdbqt +  " -U 'nphs' > out_rw_" + str(n) + ".log"
        cmd_pw ="$VINADIR/vina --receptor " + propdbqt + " --ligand " + wpdbqt + "  --score_only --log score_RW_" + str(n) + ".txt > out_pw" + str(n) + ".log"
        try:
            os.system(cmd2)
            os.system(cmd_pw)
        except:
            os.system("touch FAIL.log")

    os.chdir(olddir)


def get_result_PW(fn,out_PW): 
    ''' Get PW result '''

    olddir = os.getcwd()
    os.chdir("vina_RW")
    num = len([file for file in os.listdir(".") if file.startswith("score_RW_") and file.endswith(".txt")])
    for i in range(num):
        try:
            ### vina in linux
            value = [float(line.split()[1]) for line in open("score_RW_" + str(i) + ".txt") if line[0:8] == "Affinity"][0]
        except:
            ### vina in os
            value = [float(line.split(":")[1]) for line in open("score_RW_" + str(i) + ".txt") if line[0:8] == "Affinity"][0]
        out_PW.write(fn + "," + str(i) + "," + str(value)  + "\n")
    out_PW.close()
    os.chdir(olddir)


def get_RW_final(fn,out,out_total):
    ''' Get RW which satisfies Vina score requirement '''

    value_PW = 0.00
    m = 0
    index = []
    for line in open("Epw.csv"):
        if line.split(",")[0] == fn:
            PW_score = float(line.split(",")[2])
            if PW_score < 0:
                index.append(line.split(",")[1])
                out.write(fn + "," + ",".join(line.split(",")[1:]))
                value_PW += float(PW_score)*-0.73349
                m = m + 1
    out.close()
    out_total.write(fn + " " + str(m) + " " + str(value_PW)  + "\n")
    
    return index

def get_waterfile(fn, pro, index):
    ''' Get final RW water file and protein file '''
    
    ### get final rw water file ###
    out = open("RW_total.pdb","w")
    for n, line in enumerate(open("water_index.txt")):
        if str(n) in index:
            line = line.rstrip()
            if "." in line:
                water_file = "RW_" + "_".join(line.rstrip().split(".")) + ".pdb"
            else:
                water_file = "RW_" + line.rstrip() + "_chain.pdb"
            lines = [ l for l in open(water_file) if (l[0:6] in ["ATOM  ","HETATM"]) and (l[0:3] != "END")]
            out.write("".join(lines))
    out.close()

    out = open("RW_total.pdb")
    count_lines = int(len(out.readlines())/3)
    out.close()

    if os.stat("RW_total.pdb").st_size == 0:
        print('No RW has been found in ' + fn)
    else:
        print(str(count_lines) + " RW have been saved in " + fn  + "_protein_RW.pdb")
    
    ### get final rw protein file ###
    in_rec = open(pro)
    out_rec = open(fn + "_protein_RW.pdb","w")
    lines = [line for line in in_rec if line[0:6] in ["ATOM  ", "HETATM"]]
    out_rec.write("".join(lines))
    in_water = open("RW_total.pdb")
    lines = [line for line in in_water]
    out_rec.write("".join(lines))
    out_rec.write("END")
    out_rec.close()


def get_Crw(fn,inprot,inwater,datadir):
    """
    Get Crw
    
    :param fn: input index
    :param inprot: inpro with only protein
    :param inwater: inpro with both protein and water
    :param datadir: directory for input 

    """

    os.chdir(datadir)
    pro = os.path.join(datadir,inprot)
    water = os.path.join(datadir,inwater)
    os.system("mkdir RW")
    os.chdir("RW")
    get_RW(fn,water)
    get_water(fn,water)
    addH(fn)
    runVina(fn,pro)
    out_PW = open("Epw.csv","w")
    get_result_PW(fn,out_PW)
    out  = open("RW.csv","w")
    out_total = open("RW_total.csv","w")
    index = get_RW_final(fn,out,out_total)
    get_waterfile(fn,pro,index)
    os.system("cp " + fn + "_protein_RW.pdb ../")
    os.chdir(datadir)
    os.system("rm -r RW*")


