#-----------------------------------------------------------------------------
# SASA Feature
#-----------------------------------------------------------------------------
import os
import sys
import fileinput
from DXGB.featureSASA import sasa

def cal_SASA(out,fn,lig,pro,datadir):
    """
    Calculate SASA features
    
    :param out: outfile object (file can be wrote)
    :param fn: input index
    :param lig: inlig file
    :param pro: inpro file 
    :param datadir: directory for input 

    return SASA 30 features for current input index will be wrote into outfile
    """
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


