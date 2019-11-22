#-----------------------------------------------------------------------------
# Combine Features
#-----------------------------------------------------------------------------
import pandas as pd
import os

def read_file(Vina58,SASA,dE,Water,Ion):
    combine = "pdb"
    f_V58 = pd.read_csv(Vina58,dtype={"pdb":str})
    f_V58["vina"] = f_V58["vina"] * -0.73349
    f_SASA = pd.read_csv(SASA,dtype={"pdb":str})
    f = pd.merge(f_V58,f_SASA, on = combine)
    if dE == None:
        f["dE_global"] = -300
        f["RMSD_global"] = 300
    else:
        f_dE = pd.read_csv(dE,dtype={"pdb":str})
        f = pd.merge(f,f_dE, on = combine)
    if Water == None:
        f["Nbw"] = 0
        f["Epw"] = 0
        f["Elw"] = 0
    else:
        f_water = pd.read_csv(Water,dtype={"pdb":str})
        f = pd.merge(f,f_water, on = combine)
    if Ion == None:
        f["Ni"] = 0
    else:
        f_ion = pd.read_csv(Ion,dtype={"pdb":str})
        f = pd.merge(f,f_ion, on = combine)
    
    f["pdb"] = f["pdb"].astype(str) 
    
    return f


def combine(datadir,data_type = ""):
    d_wat = {"":None, "_BW": os.path.join(datadir,"Feature_BW_BW.csv"), "_RW":os.path.join(datadir,"Feature_BW_RW.csv"), "_min":None, 
            "_min_RW":os.path.join(datadir,"Feature_BW_min_RW.csv"), 
            "_min_BW":os.path.join(datadir,"Feature_BW_min_BW.csv"), 
            "_min_PW":os.path.join(datadir,"Feature_BW_min_PW.csv")}
    Vina58 = os.path.join(datadir,"Vina58" + data_type + ".csv")
    SASA = os.path.join(datadir,"SASA" + data_type + ".csv")
    dE = os.path.join(datadir,"dE_RMSD.csv")
    Water = d_wat[data_type]
    Ion = os.path.join(datadir, "Num_Ions" + data_type + ".csv")
    outfile = "Input" + data_type + ".csv"
    f = read_file(Vina58,SASA,dE,Water,Ion)
    f.to_csv(os.path.join(datadir,outfile),index = False)
    