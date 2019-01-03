
import os
import get_pdbinfo
import numpy as np
import sys
from software_path import path_obabel

obable = path_obabel()

def get_Ions(fn,lig,pro,infile):
    outfile = open(infile,"w")

    lig = get_pdbinfo.pdbinfo(fn,file = lig)
    lig_atoms = lig.getPolarAtoms()

    pro = get_pdbinfo.pdbinfo(fn,file = pro)
    pro_ions = pro.getIons()

    ion_coord = get_pdbinfo.pdbinfo(fn,lines = pro_ions).getCoords()
    lig_coord = get_pdbinfo.pdbinfo(fn,lines = lig_atoms).getCoords()

    ion_coord = np.expand_dims(ion_coord, 1)
    lig_coord = np.expand_dims(lig_coord,0)

    distance = np.linalg.norm(ion_coord - lig_coord, axis = 2)
    distance_min = np.min(distance, axis = 1)

    ion_index = []
    for idx, i in enumerate(distance_min):
        if i < 3.5:
            ion_index.append(idx)

    for i in ion_index:
        ion_line = pro_ions[i]
        ion_distance = distance[i]
        ion_chain = get_pdbinfo.chid(ion_line)
        ion_name = get_pdbinfo.atmn(ion_line).strip()

        ion_chain = get_pdbinfo.chid(ion_line)
        if ion_chain != " ":
            ion_idx = str(int(get_pdbinfo.resi(ion_line))) + "." + ion_chain 
        else:
            ion_idx = str(int(get_pdbinfo.resi(ion_line)))
            
        for idx, d in enumerate(ion_distance):
            if d < 3.5:
                lig_line = lig_atoms[idx]
                lig_name = get_pdbinfo.atmn(lig_line).strip()
                outline = fn + "," + ion_idx + "," + ion_name + "," + lig_name + "," + str(round(d,2)) + "\n"
                outfile.write(outline)
    outfile.close()

    return None

def get_num(fn, infile):
    ion = []
    for line in open(infile):
        ion.append(line.split(",")[1] + line.split(",")[2])
    ionset = set(ion)
    return len(ionset)

def cal_Ni(outfile,fn, inprot, inlig, datadir):
    pro = os.path.join(datadir, inprot)
    lig = os.path.join(datadir,inlig)
    infile = os.path.join(datadir,"Ion_infor.dat")
    get_Ions(fn,lig,pro,infile)
    num = get_num(fn,infile)
    #os.system("rm " + infile)
    outfile.write(fn + "," + str(num) + "\n")

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
    outfile = open(datadir + "Num_Ions.csv","w")
    outfile_RW = open(datadir + "Num_Ions_RW.csv","w")
    outfile_min = open(datadir + "Num_Ions_min.csv","w")
    outfile_min_RW = open(datadir + "Num_Ions_min_RW.csv","w")
    outfile.write("pdb,Ni\n")
    outfile_RW.write("pdb,Ni\n")
    outfile_min.write("pdb,Ni\n")
    outfile_min_RW.write("pdb,Ni\n")
    for fn in pdblist:
        inpro = fn + "_protein.pdb"
        inpro_RW = fn + "_protein_RW.pdb"
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
            os.system(obable + " -ipdbqt " + datadir + inlig + " -omol2 -O " + datadir+ inlig_out)
            inlig = inlig_out

        else:
            print("wrong ligand input file format, it should be mol2, sdf, or pdbqt")
        cal_Ni(outfile,fn,inpro,inlig,datadir)
        cal_Ni(outfile_RW,fn,inpro_RW,inlig,datadir)
        inlig_min = fn + "_lig_min.pdb"
        cal_Ni(outfile_min,fn,inpro,inlig_min,datadir)
        inlig_min_RW = fn + "_lig_min_RW.pdb"
        cal_Ni(outfile_min_RW,fn,inpro_RW,inlig_min_RW,datadir)

    outfile.close()
    outfile_RW.close()
    outfile_min.close()
    outfile_min_RW.close()

    return None

if __name__ == "__main__":
    #main()
    outfile = open("/Users/jianinglu1/Documents/script/deltaXGB_linux/Feature/test/Num_Ions.csv","w")
    fn = "4q90"
    inprot = fn + "_protein.pdb"
    inlig = fn + "_ligand.pdb"
    datadir = "/Users/jianinglu1/Documents/script/deltaXGB_linux/Feature/test"
    cal_Ni(outfile,fn, inprot, inlig, datadir)

