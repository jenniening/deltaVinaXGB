import os
import pandas as pd
from software_path import path_R

R = path_R()
def convert_RF20(datadir,infile,outfile):
    print(datadir + infile)
    infile = pd.read_csv(datadir + infile)
    features_vina = ["vina1", "vina3","vina53","vina55","vina54","vina56","vina4","vina52","vina58","vina48"]
    features_SASA = ["P2.P","P2.N","P2.DA","P2.D","P2.A","P2.AR","P2.H","P2.PL","P2.HA","P2.SA"]
    features_vina_change = ["vina1","vina3","vina4","vina52","vina48"]
    features_name = ["F1","F2","F3","F4","F5","F6","F7","F8","F9","F10","F11","F12","F13","F14","F15","F16","F17","F18","F19","F20"]
    features_list = ["pdb","vina"] + features_vina + features_SASA
    newfile = infile[features_list]
    newfile[features_vina_change] = newfile[features_vina_change] * -0.73349
    newfile.columns = ["pdb","vina"] + features_name
    newfile.to_csv(datadir + outfile,index = False)

    return None

def get_RF20(datadir, infile, outfile):
    filename = "get_RF20_new.R"
    newfile = open(filename,"w")
    lines = [line for line in open("get_RF20.R")]
    newfile.write("".join(lines[0:6]))
    newfile.write("infn = '" + datadir + infile + "'\n")
    newfile.write("outfn = '" + datadir + outfile + "'\n")
    newfile.write("".join(lines[11:]))
    newfile.close()

    os.system(R + " CMD BATCH " + filename)

    return None

def RF20_main(datadir,infile, outfile):
    convert_RF20(datadir,infile,"RF_input.csv")
    get_RF20(datadir,"RF_input.csv",outfile)

    return None



