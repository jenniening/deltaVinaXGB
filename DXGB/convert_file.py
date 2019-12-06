#-----------------------------------------------------------------------------
# Convert File for deltaVinaRF20
#-----------------------------------------------------------------------------
import sys
import os
import pandas as pd
from DXGB.utils import get_tool

project = get_tool("project")
Rscript = get_tool("Rscript")

def convert_RF20(infile,outfile):
    """
    Change several features into pKd value, since original deltaVinaRF20 trained using converted value 

    :param infile: initial feature file
    :param outfile: output file name with converted features

    """

    infile = pd.read_csv(infile,dtype = {"pdb":str})
    features_vina = ["vina1", "vina3","vina53","vina55","vina54","vina56","vina4","vina52","vina58","vina48"]
    features_SASA = ["P2.P","P2.N","P2.DA","P2.D","P2.A","P2.AR","P2.H","P2.PL","P2.HA","P2.SA"]
    features_vina_change = ["vina1","vina3","vina4","vina52","vina48"]
    features_list = ["pdb","vina"] + features_vina + features_SASA
    newfile = infile[features_list]
    newfile[features_vina_change] = newfile[features_vina_change] * -0.73349
    newfile.to_csv(outfile,index = False)

def get_RF20(RF20da, infile, outfile):
    """
    Get deltaRF20 score 
    
    :param RF20da: deltaVinaRF20 model 
    :param infile: input feature file
    :param outfile: output score file 

    """
    RFfile = open("get_RF20.R","w")
    
    RFfile.write("library(randomForest)\n" + "load('" + RF20da + "')\n" + 
                 "args <- commandArgs(trailingOnly = TRUE)\ninfn = args[1]\noutfn = args[2]\n" + 
                 "df = read.table(infn, header=T, stringsAsFactors = F, sep=',')\n" + 
                 "feats = df[3:22]\n" + "pred = predict(rffit, newdata = feats) + df$vina\n" +
                 "output = data.frame(pdb = df$pdb, RF20 = pred)\n" + 
                 "write.table(output, outfn, sep=',', row.names = F, quote = F)")
    RFfile.close()
    os.system(Rscript + " get_RF20.R " + infile + " " + outfile)

	
def RF20_main(datadir,infile, outfile, RFmodel="RF20_rm2016.rda"):
    """[summary]
    
    :param datadir: [description]
    :type datadir: [type]
    :param infile: [description]
    :type infile: [type]
    :param outfile: [description]
    :type outfile: [type]
    """
    
    olddir = os.getcwd()
    os.chdir(datadir)
    if RFmodel not in os.listdir(datadir):
        os.system("cp " + os.path.join(project, RFmodel) + " .")
    RF20da = RFmodel
    convert_RF20(infile,"RF_input.csv")
    get_RF20(RF20da, "RF_input.csv",outfile)
    os.chdir(olddir)
