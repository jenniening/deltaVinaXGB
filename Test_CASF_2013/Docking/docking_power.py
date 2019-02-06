import numpy as np
import pandas as pd

def dockPower(dk13):
    """
    dk13 columns: pdb, Index, RMSD, sc,

    """

    dk13['rank'] = dk13.groupby(['pdb'])['sc'].rank(method = 'min', ascending = False)
    dk13.head()

    sr = []
    for i in [1]:
        # count the cases rank < (2,3,4) and RMSD < 2.0
        # total 195 proteins
        sr.append(str(len(dk13.loc[(dk13['rank'] == i) & (dk13['RMSD'] < 2.0)]['pdb'].unique())))
        print(dk13.loc[(dk13['rank'] == i) & (dk13['RMSD'] > 2.0)]['pdb'].unique())
    return sr


def get_Docking(infile):
    """
    Get the docking success rate of model
    :param dk13: CASF-2013 docking test data
    :param model: model_index
    :return dc: list of docking success rates for the model
    """
    score = pd.read_csv(infile)
    
    dk_vina = dk13["vina"]
        dk_pred = dk_total.sum(axis = 1)/1  + dk_vina
        dk_new = pd.concat([dk13[["pdb","index","RMSD"]],dk_pred],axis = 1)
        dk_new.columns = ['pdb','Index','RMSD','sc']
        dc.append(dockPower(dk_new)[0])

    return dc

