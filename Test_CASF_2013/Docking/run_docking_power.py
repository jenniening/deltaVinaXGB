import os

def run_docking(score_dir, out_dir, relation):
    os.system("python docking_power.py " + score_dir +  " " + out_dir + " " + score_dir.split("/")[-1] + "_docking.out " + relation)

if __name__ == "__main__":
    datadir = "/scratch/jl7003/deltaVinaXGB_develop/Test_CASF_2013/Docking/"
    #datadir = "/beegfs/jl7003/CASF-2016/power_scoring/other_score"
    #model_name = "model_17"
    #model_name = "Vina"
    relation = "positive"
    for n in ["34"]:
        for i in ["C","Co"]:
            score_dir = datadir + "our_score_" + i + "_final_model" + n  + "_nocry"
            #score_dir = datadir + "Vina_" + i + "_final_nocry"
            out_dir = datadir + "result/"
            run_docking(score_dir,out_dir,relation)
