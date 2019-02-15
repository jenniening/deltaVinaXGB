import os

def run_ranking(score_name, score_dir, relation):
    os.system("python ranking_power.py " + score_dir +  " " + score_name + " " + score_name.split(".")[0] + "_ranking.out " + relation)

if __name__ == "__main__":
    datadir = "/scratch/jl7003/deltaVinaXGB_develop/Test_CASF_2013/Ranking/our_score_0129"
    #datadir = "/beegfs/jl7003/CASF-2016/power_scoring/other_score"
    model_name = "model_37"
    #model_name = "RF_remove2016"
    relation = "negative"
    for i in ["C","Co","Crwo","Cbwo"]:
        score_dir = datadir + "/" + model_name + "/"
        score_name = "CASF-2013_Scoring_score_" + i + ".csv"
        run_ranking(score_name, score_dir, relation)
