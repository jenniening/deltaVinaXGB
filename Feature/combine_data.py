import pandas as pd
import os

def read_file(Vina58,SASA,dE,Water,Ion,Core,Side,NumFrags,decoy):
    if decoy:
        combine = ["pdb","idx"]
    else:
        combine = "pdb"
    f_V58 = pd.read_csv(Vina58,dtype={"pdb":str})
    f_V58["vina"] = f_V58["vina"] * -0.73349
    f_SASA = pd.read_csv(SASA,dtype={"pdb":str})
    f = pd.merge(f_V58,f_SASA, on = combine)
    if dE == None:
        f["dE_global"] = -300
        f["RMSD_global"] = 300
        f["number0"] = 300
        f["number1"] = 300
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
    f_core = pd.read_csv(Core,dtype={"pdb":str})
    f_side = pd.read_csv(Side,dtype={"pdb":str})
    f_numfrag = pd.read_csv(NumFrags,dtype={"pdb":str})

    f_total = pd.merge(f_core,f_side,on = combine)
    f = pd.merge(f, f_total, on = combine)
    f = pd.merge(f,f_numfrag, on = combine)
    f["pdb"] = f["pdb"].astype(str) 
    if decoy:
        f["idx"] = f["idx"].astype(str)
        assert f.shape[1] == 218, "Feature Calculation Failed"
    else:
        assert f.shape[1] == 217, "Feature Calculation Failed"
    
    return f


def combine(datadir,data_type = "", decoy = False):
    if not decoy:
        d_frag = {"":0,"_min":1,"_min_RW":2}
        d_wat = {"":None, "_min":None, "_min_RW":os.path.join(datadir,"Feature_BW_min_RW.csv")}
        Vina58 = os.path.join(datadir,"Vina58" + data_type + ".csv")
        SASA = os.path.join(datadir,"SASA" + data_type + ".csv")
        dE = os.path.join(datadir,"dE_RMSD.csv")
        Water = d_wat[data_type]
        Ion = os.path.join(datadir, "Num_Ions" + data_type + ".csv")
        NumFrags = os.path.join(datadir,"NumFrags.csv")
        Core = os.path.join(datadir,"Frags/Vina_core_" + str(d_frag[data_type]) + ".csv")
        Side = os.path.join(datadir,"Frags/Vina_side_" + str(d_frag[data_type]) + ".csv")
        outfile = "Input" + data_type + ".csv"
        f = read_file(Vina58,SASA,dE,Water,Ion,Core,Side,NumFrags,decoy)
        f.to_csv(os.path.join(datadir,outfile),index = False)
    else:
        d_wat = {"":None, "_RW":os.path.join(datadir,"Feature_BW_decoys" + data_type + ".csv")}
        Vina58 = os.path.join(datadir,"Vina58_decoys" + data_type + ".csv")
        SASA = os.path.join(datadir,"SASA_decoys" + data_type + ".csv")
        dE = os.path.join(datadir,"dE_RMSD_decoys" + data_type + ".csv")
        Water = d_wat[data_type]
        Ion = os.path.join(datadir, "Num_Ions_decoys" + data_type + ".csv")
        NumFrags = os.path.join(datadir,"NumFrags_decoys" + data_type + ".csv")
        Core = os.path.join(datadir,"Vina_core_decoys" + data_type + ".csv")
        Side = os.path.join(datadir,"Vina_side_decoys" + data_type + ".csv")
        outfile = "Input_decoys" + data_type + ".csv"
        f = read_file(Vina58,SASA,dE,Water,Ion,Core,Side,NumFrags,decoy)
        f.to_csv(os.path.join(datadir,outfile),index = False)



        
    



if __name__ == "__main__":

    testdir = '/Users/jianinglu1/Documents/script/deltavinaXGB/Test/'
    combine(testdir)
