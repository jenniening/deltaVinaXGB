import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
import sys
import os


def get_Ranking(score, relation = "positive"):
    """
    Get the ranking performance of model
    :param model: model_index
    :return High_Rank_list, Low_Rank_list: Correct order of 3 ligands and only choose top1 correctly
    """
    High_Rank_list = []
    Low_Rank_list = []
    score = pd.read_csv(score, sep = " ")
    if relation == "negative":
        score["score"] = -score["score"]

    """ Ranking Performance """
    refdata = pd.read_csv("Score_Rank.dat",sep = " ")
    refdata.columns = ["#code", "Prot", "lig", "pKd"]
    newdf = pd.merge(refdata,score,on = "#code")
    newdf["rank_org"] = newdf.groupby("Prot")["pKd"].rank(method = "min",ascending = False)
    newdf["rank_pre"] = newdf.groupby("Prot")["score"].rank(method = "min", ascending = False)
    ### High ###
    data = newdf[newdf["rank_org"] == newdf["rank_pre"]].groupby("Prot").size()
    success = len([ i for i in data if i == 3])
    print("High Rank:" + str(success) + "," + str(round(success/65*100,1)))
    High_Rank = str(round(success/65 *100,1))
    High_Rank_list.append(High_Rank)
    ### Low ###
    data = newdf[newdf["rank_org"] == newdf["rank_pre"]][newdf["rank_org"].isin([1])]
    print("Low Rank:" + str(data.shape[0]) + "," + str(round(data.shape[0]/65 *100,1)))
    Low_Rank = str(round(data.shape[0]/65 *100,1))
    Low_Rank_list.append(Low_Rank)
    return High_Rank,Low_Rank

if __name__ == "__main__":  
    args = sys.argv[1:]
    datadir = args[0]
    score = args[1]
    outfile = open(os.path.join(datadir,args[2]), "w")
    outfile.write("HR, LR\n")
    relation = args[3]
    HR, LR = get_Ranking(os.path.join(datadir,score), relation = relation)
    outfile.write(str(HR) + "," + str(LR))
    outfile.close()
    #score = pd.read_csv(score, sep = " ")
    #df = pd.merge(exp, score, on = "#code")

