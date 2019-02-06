import os

def run_scoring(score_name, score_dir):
    os.system("python scoring_power.py " + score_dir +  " " + score_name + " " + score_name + ".out")

if __name__ == "__main__":
    datadir = "/beegfs/jl7003/CASF-2016/power_scoring/our_score_0129"
    #datadir = "/beegfs/jl7003/CASF-2016/power_scoring/other_score"
    model_name = "model_16"
    #model_name = "RF"
    for i in ["C","Co","Crwo","Cbwo"]:
        score_dir = datadir + "/" + model_name + "/CASF-2013_Scoring_score_" + i + ".csv"
        score_name = datadir + "/" + model_name + "/" + score_dir.split("/")[-1].split(".")[0]
        run_scoring(score_name, score_dir)
