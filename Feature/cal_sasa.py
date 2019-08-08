#-----------------------------------------------------------------------------
# SASA Feature
#-----------------------------------------------------------------------------
import os
import sys
import fileinput
from Feature.featureSASA import sasa

def cal_SASA(out,fn,lig,pro,datadir):
    # fn: pdbid
    # lig: ligand part
    # pro: protein part
    pro = os.path.join(datadir,pro)
    lig = os.path.join(datadir,lig)
    sasa_features = sasa(datadir,pro,lig)
    sasa_com = sasa_features.sasa
    sasa_pro = sasa_features.sasa_pro
    sasa_lig = sasa_features.sasa_lig

    out.write(fn + "," +  ",".join([str(round(i,2) )for i in sasa_com]) + "," + ",".join([str(round(i,2) )for i in sasa_lig]) + "," + ",".join([str(round(i,2) )for i in sasa_pro]) + "\n")
    return None

if __name__ == "__main__":
    datadir  =  "/Users/jianinglu1/Documents/GitHub/deltaVinaXGB_develop/Feature/test"
    fn = "3c2f"
    lig = fn + "_ligand.pdb"
    pro = fn + "_protein.pdb"
    out = open(datadir + "/out_sasa.csv","w")
    cal_SASA(out,fn,lig,pro,datadir)
    out.close()


