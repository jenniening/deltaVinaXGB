import os

def run_docking(score_dir, out_dir, relation):
    os.system("python docking_power.py " + score_dir +  " " + out_dir + " " + score_dir.split("/")[-1] + "_docking.out " + relation)

if __name__ == "__main__":
    datadir = "/scratch/jl7003/deltaVinaXGB_develop/Test_CASF_2013/Docking/"
    #datadir = "/beegfs/jl7003/CASF-2016/power_scoring/other_score"
    model_name = "model_17"
    #model_name = "RF"
    relation = "negative"
    for i in ["C","Co","Crwo","Cbwo"]:
        score_dir = datadir + "/" + model_name + "/"
        outdir_dir = datadir + "/" + model_name + "/"
        run_ranking(score_name, score_dir, relation)