import pandas as pd
import os

def read_file(Vina58,SASA,dE,Water,Ion,Core,Side,NumFrags):
    f_V58 = pd.read_csv(Vina58)
    f_V58["vina"] = f_V58["vina"] * -0.73349
    f_SASA = pd.read_csv(SASA)
    f = pd.merge(f_V58,f_SASA, on = ["pdb"])
    if dE == "none":
        f["dE_lowest"] = 0
        f["number"] = 0
    else:
        f_dE = pd.read_csv(dE)
        f = pd.merge(f,f_dE, on = ["pdb"])
    if Water == "none":
        f["Nbw"] = 0
        f["Epw"] = 0
        f["Elw"] = 0
    else:
        f_water = pd.read_csv(Water)
        f = pd.merge(f,f_water, on = ["pdb"])
    if Ion == "none":
        f["Ni"] = 0
    else:
        f_ion = pd.read_csv(Ion)
        f = pd.merge(f,f_ion, on = ["pdb"])
    f_core = pd.read_csv(Core)
    f_side = pd.read_csv(Side)
    f_numfrag = pd.read_csv(NumFrags)

    f_total = pd.merge(f_core,f_side,on = "pdb")
    f = pd.merge(f, f_total, on = "pdb")
    f = pd.merge(f,f_numfrag, on = "pdb")

    assert f.shape[1] == 217, "Feature Calculation Failed"
    
    return f


def combine(datadir,data_type = "", decoy = False):
    if not decoy:
        d_frag = {"":0,"_min":1,"_min_RW":2}
        d_wat = {"":none, "_min":none, "_min_RW":os.path.join(datadir,"Feature_BW_min_RW.csv")}
        if data_type == "":
            Vina58 = os.path.join(datadir,"Vina58" + data_type + ".csv")
            SASA = os.path.join(datadir,"SASA" + data_type + ".csv")
            dE = os.path.join(datadir,"dE_RMSD.csv")
            Water = d_wat[data_type]
            Ion = os.path.join(datadir, "Num_Ions" + data_type + ".csv")
            NumFrags = os.path.join(datadir,"NumFrags.csv")
            Core = os.path.join(datadir,"Frags/Vina_core_" + str(d_frag[data_type] + ".csv"))
            Side = os.path.join(datadir,"Frags/Vina_side_" + str(d_frag[data_type] + ".csv"))
            outfile = "Input" + data_type + ".csv"
            f = read_file(Vina58,SASA,dE,Water,Ion,Core,Side,NumFrags)
            f.to_csv(datadir + outfile,index = False)
    



if __name__ == "__main__":

    testdir = '/Users/jianinglu1/Documents/script/deltavinaXGB/Test/'
    combine(testdir)
