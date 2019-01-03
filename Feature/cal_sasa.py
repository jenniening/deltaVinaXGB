"""
SASA 30 Features
"""
__author__ = "Jianing Lu"
__copyright__ = "Copyright 2018, NYU"
__license__ = ""

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import featureSASA
from featureSASA import sasa
import fileinput
import os
import sys
from software_path import path_obabel

obable = path_obabel()

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------
def cal_SASA(out,fn,lig,pro,datadir):
    ## fn: pdbid
    ## lig: ligand part
    ## pro: protein part
    pro = os.path.join(datadir,pro)
    lig = os.path.join(datadir,lig)
    sasa_features = sasa(datadir,pro,lig)
    sasa_com = sasa_features.sasa
    sasa_pro = sasa_features.sasa_pro
    sasa_lig = sasa_features.sasa_lig

    out.write(fn + "," +  ",".join([str(round(i,2) )for i in sasa_com]) + "," + ",".join([str(round(i,2) )for i in sasa_lig]) + "," + ",".join([str(round(i,2) )for i in sasa_pro]) + "\n")

    return None



def main():
    args = sys.argv[1:]
    if args[-1] == "file":
        pdbfile = open('%s'%(sys.argv[1] + sys.argv[2]),'r')
        pdblist = []
        for i in pdbfile.readlines():
            pdblist.append(i[0:4])
    else:
        pdblist = []
        pdblist.append(sys.argv[2])
    datadir = sys.argv[1]
    out_SASA = open(datadir + "SASA.csv","w")
    out_SASA_RW = open(datadir + "SASA_RW.csv","w")
    out_SASA_min = open(datadir + "SASA_min.csv","w")
    out_SASA_min_RW = open(datadir + "SASA_min_RW.csv","w")
    f_type = ["P","N","DA","D","A","AR","H","PL","HA","SA"]
    f_com = list("P2." + i for i in f_type)
    f_dl = list("P2dl." + i for i in f_type)
    f_dp = list("P2dp." + i for i in f_type)
    out_SASA.write("pdb," + ",".join(f_com) + "," + ",".join(f_dl) + "," + ",".join(f_dp) + "\n")
    out_SASA_RW.write("pdb," + ",".join(f_com) + "," + ",".join(f_dl) + "," + ",".join(f_dp) + "\n")
    out_SASA_min.write("pdb," + ",".join(f_com) + "," + ",".join(f_dl) + "," + ",".join(f_dp) + "\n") 
    out_SASA_min_RW.write("pdb," + ",".join(f_com) + "," + ",".join(f_dl) + "," + ",".join(f_dp) + "\n")
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
            os.system(obable + " -ipdbqt " + datadir + inlig + " -omol2 -O " + datadir + inlig_out)
            inlig = inlig_out

        else:
            print("wrong ligand input file format, it should be mol2, sdf, or pdbqt")
        cal_SASA(out_SASA,fn,inlig,inpro,datadir)
        cal_SASA(out_SASA_RW,fn,inlig,inpro_RW,datadir)
        inlig_min = fn + "_lig_min.pdb"
        cal_SASA(out_SASA_min,fn,inlig_min,inpro,datadir)
        inlig_min_RW = fn + "_lig_min_RW.pdb"
        cal_SASA(out_SASA_min_RW,fn,inlig_min_RW,inpro_RW,datadir)
    out_SASA.close()
    out_SASA_RW.close()
    out_SASA_min.close()
    out_SASA_min_RW.close()

    return None

if __name__ == "__main__":

    #main()
    datadir  =  "/Users/jianinglu1/Documents/GitHub/deltaVinaXGB_develop/Feature/test"
    fn = "3c2f"
    lig = fn + "_ligand.pdb"
    pro = fn + "_protein.pdb"
    out = open(datadir + "/out_sasa.csv","w")
    cal_SASA(out,fn,lig,pro,datadir)
    out.close()


