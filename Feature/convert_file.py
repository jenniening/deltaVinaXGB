#-----------------------------------------------------------------------------
# Convert File for deltaVinaRF20
#-----------------------------------------------------------------------------
import sys
import os
import pandas as pd

if sys.platform == "linux":
    from Feature.software_path_linux import *
elif sys.platform == "darwin":
    from Feature.software_path_mac import *
	
Rscript = path_Rscript()

def convert_RF20(infile,outfile,decoy):
    if decoy:
        infile = pd.read_csv(infile,dtype = {"pdb":str, "idx":str})
    else:
        infile = pd.read_csv(infile,dtype = {"pdb":str})
    features_vina = ["vina1", "vina3","vina53","vina55","vina54","vina56","vina4","vina52","vina58","vina48"]
    features_SASA = ["P2.P","P2.N","P2.DA","P2.D","P2.A","P2.AR","P2.H","P2.PL","P2.HA","P2.SA"]
    features_vina_change = ["vina1","vina3","vina4","vina52","vina48"]
    features_name = ["F1","F2","F3","F4","F5","F6","F7","F8","F9","F10","F11","F12","F13","F14","F15","F16","F17","F18","F19","F20"]
    if decoy:
        features_list = ["pdb","idx","vina"] + features_vina + features_SASA
    else:
        features_list = ["pdb","vina"] + features_vina + features_SASA
    newfile = infile[features_list]
    newfile[features_vina_change] = newfile[features_vina_change] * -0.73349
    if decoy:
        newfile.columns = ["pdb","idx","vina"] + features_name
    else:
        newfile.columns = ["pdb","vina"] + features_name
    newfile.to_csv(outfile,index = False)
    return None
	
def get_RF20(infile, outfile, decoy):
    RFfile = open("get_RF20.R","w")
    RF20da = path_RF20da()
    if decoy:
        RFfile.write("library(randomForest)\n" + "load('" + RF20da + "')\n" +
                 "args <- commandArgs(trailingOnly = TRUE)\ninfn = args[1]\noutfn = args[2]\n" + 
                 "df = read.table(infn, header=T, stringsAsFactors = F, sep=',')\n" +
                 "feats = df[4:23]\n" + "pred = predict(rffit, newdata = feats) + df$vina\n" +
                 "output = data.frame(pdb = df$pdb, idx = df$idx, RF20 = pred)\n" +
                 "write.table(output, outfn, sep=',', row.names = F, quote = F)")
    else:
        RFfile.write("library(randomForest)\n" + "load('" + RF20da + "')\n" + 
                 "args <- commandArgs(trailingOnly = TRUE)\ninfn = args[1]\noutfn = args[2]\n" + 
                 "df = read.table(infn, header=T, stringsAsFactors = F, sep=',')\n" + 
                 "feats = df[3:22]\n" + "pred = predict(rffit, newdata = feats) + df$vina\n" +
                 "output = data.frame(pdb = df$pdb, RF20 = pred)\n" + 
                 "write.table(output, outfn, sep=',', row.names = F, quote = F)")
    RFfile.close()
    os.system(Rscript + " get_RF20.R " + infile + " " + outfile)
    return None
	
def RF20_main(datadir,infile, outfile, decoy):
    olddir = os.getcwd()
    os.chdir(datadir)
    convert_RF20(infile,"RF_input.csv", decoy)
    get_RF20("RF_input.csv",outfile, decoy)
    os.chdir(olddir)
