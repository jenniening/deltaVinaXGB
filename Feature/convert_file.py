import sys
import os
import pandas as pd
if sys.platform == "linux":
    from software_path_linux import path_Rscript
    from software_path_linux import path_RF20da
elif sys.platform == "darwin":
    from software_path_mac import path_Rscript
    from software_path_mac import path_RF20da
	
Rscript = path_Rscript()
def convert_RF20(infile,outfile):
    infile = pd.read_csv(infile)
    features_vina = ["vina1", "vina3","vina53","vina55","vina54","vina56","vina4","vina52","vina58","vina48"]
    features_SASA = ["P2.P","P2.N","P2.DA","P2.D","P2.A","P2.AR","P2.H","P2.PL","P2.HA","P2.SA"]
    features_vina_change = ["vina1","vina3","vina4","vina52","vina48"]
    features_name = ["F1","F2","F3","F4","F5","F6","F7","F8","F9","F10","F11","F12","F13","F14","F15","F16","F17","F18","F19","F20"]
    features_list = ["pdb","vina"] + features_vina + features_SASA
    newfile = infile[features_list]
    newfile[features_vina_change] = newfile[features_vina_change] * -0.73349
    newfile.columns = ["pdb","vina"] + features_name
    newfile.to_csv(outfile,index = False)
    return None
	
def get_RF20(infile, outfile):
    RFfile = open("get_RF20.R","w")
    RF20da = path_RF20da()
    RFfile.write("library(randomForest)\n" + "load('" + RF20da + "')\n" + 
                 "args <- commandArgs(trailingOnly = TRUE)\ninfn = args[1]\noutfn = args[2]\n" + 
                 "df = read.table(infn, header=T, stringsAsFactors = F, sep=',')\n" + 
                 "feats = df[3:22]\n" + "pred = predict(rffit, newdata = feats) + df$vina\n" +
                 "output = data.frame(pdb = df$pdb, RF20 = pred)\n" + 
                 "write.table(output, outfn, sep=',', row.names = F, quote = F)")
    RFfile.close()
    os.system(Rscript + " get_RF20.R " + infile + " " + outfile)
    return None
	
def RF20_main(datadir,infile, outfile):
    olddir = os.getcwd()
    os.chdir(datadir)
    convert_RF20(infile,"RF_input.csv")
    get_RF20("RF_input.csv",outfile)
    os.chdir(olddir)
