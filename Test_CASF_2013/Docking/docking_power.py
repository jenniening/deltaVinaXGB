import numpy as np
import pandas as pd
import os
import sys

def get_Docking(pdblist, relation = "positive"):
    """
    dk13 columns: #code, score

    """
    sr_1 = 0
    sr_2 = 0
    sr_3 = 0 
    correct_list = []
    for fn in pdblist:
        dk13 = pd.read_csv(fn + "_score.dat", sep = " ")
        if relation = "negative":
            dk13["score"] = -dk13["score"]
        rmsd = pd.read_csv(fn + "_rmsd.dat", sep = " ", header = None)
        rmsd.columns = ["#code","rmsd_1","rmsd"]
        dk13 = pd.merge(dk13,rmsd[["#code","rmsd"]], on = "#code")

        dk13['rank'] = dk13.groupby(['#code'])['score'].rank(method = 'min', ascending = False)
        # count the cases rank < (2,3,4) and RMSD < 2.0
        # total 195 proteins
        if len(dk13.loc[(dk13['rank'] < 2) & (dk13['RMSD'] < 2.0)]['#code'].unique()) == 1:
            sr_1 += 1
        if len(dk13.loc[(dk13['rank'] < 3) & (dk13['RMSD'] < 2.0)]['#code'].unique()) == 1:
            sr_2 += 1
        if len(dk13.loc[(dk13['rank'] < 4) & (dk13['RMSD'] < 2.0)]['#code'].unique()) == 1:
            sr_3 += 1
        correct_list.append(dk13.loc[(dk13['rank'] < 1) & (dk13['RMSD'] > 2.0)]['pdb'].unique())
    sr_1 = round(sr_1 / 195 * 100, 2)
    sr_2 = round(sr_2 / 195 * 100, 2)
    sr_3 = round(sr_3 / 195 * 100, 2)

    return sr_1, sr_2, sr_3, correct_list

if __name__ == "__main__":
    args = sys.argv[1:]
    datadir = args[0]
    outdir = args[1]
    olddir = os.getcwd()
    os.chdir(datadir)
    pdblist = [fn for fn in os.listdir(datadir) if fn.endswith("_score.dat")]
    if len(pdblist) != 195:
        sys.exit("Error: input file should be prepared correctly. ")
    outfile = open(os.path.join(outdir,args[2]), "w")
    outfile.write("SR_1,SR_2,SR_3\n")
    relation = args[3]
    SR_1,SR_2,SR_3,correct_list = get_Docking(pdblist, relation = relation)
    outfile.write(str(SR_1) + "," + srt(SR_2) + "," + str(SR_3) + "\n")
    outfile.close()
    os.chdir(olddir)

