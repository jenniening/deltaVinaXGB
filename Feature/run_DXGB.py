import os
import sys
import click
import pandas as pd

if sys.platform == "linux":
    from software_path_linux import path_python
    from software_path_linux import path_obabel
elif sys.platform == "darwin":
    from software_path_mac import path_python
    from software_path_mac import path_obabel
from convert_file import RF20_main
import run_features
from run_features import run_features
import run_models
from run_models import run_model
from run_models import get_output
import click




@click.command()
@click.option("--model", default = "model_allfeatures", show_default = True, help = "model name")
@click.option("--modeldir",default = "../Model", show_default = True, help = "absolute model directory")
@click.option("--datadir", default = "../Test_linux", show_default= True, help = "absolute data directory")
@click.option("--pdbid", default = "01", show_default = True, help = "pdbid, ligand input should be pdbid_ligand.mol2 or sdf,\nprotein input should be pdbid_protein_all.pdb")
@click.option("--outfile", default = "score.csv",show_default = True, help = "output filename")
@click.option("--runfeatures",is_flag = True, show_default = True, help = "run features calculation")
@click.option("--water", default = "rw", show_default = True, help = "water type")
@click.option("--opt", default = "wo", show_default = True, help = "opt type")
@click.option("--decoy", is_flag = True, help = "decoy flag, if True, water = 'n', opt = 'n'")
@click.option("--rewrite", is_flag = True, help = "rewrite protein_RW, ligand_opt, generated confs or not")
@click.option("--average",is_flag = True, help = "average for 10 models")
@click.option("--modelidx", default = "1", show_default = True, help = "model index")

def main(model, modeldir, datadir, pdbid, outfile, runfeatures, water, opt, decoy, rewrite, average, modelidx):
    datadir = os.path.realpath(datadir)
    print("pdb index: " + pdbid  )
    print("file directory: " + datadir)
    print("output filename : " + outfile)
    olddir = os.getcwd()
    if runfeatures:
        run_features(datadir,pdbid, water_type = water, opt_type = opt, ligand_type = "wo", rewrite = rewrite, decoy = decoy)
        os.chdir(olddir)



    if decoy:
        opt = "n"

    if opt == "wo":
        data_type = ["","_min","_min_RW"]
    elif opt == "o":
        data_type = ["","_min"]
    else:
        data_type = [""]

    out = []
    if decoy:
        datadir = os.path.join(datadir,"decoy")
    for i in data_type:
        inf = "Input" + i + ".csv"
        test_new = run_model(inf,datadir,i,model_dir = modeldir, model_name = model, average = average, model_index = modelidx)
        outRF = "RF" + i + ".csv"
        RF20_main(datadir,inf,outRF)
        outRF_new = open(os.path.join(datadir,"RF" + i + "_new.csv"),"w")
        outRF_new.write("pdb,RF20" + i + "\n")
       	lines = [line for line in open(os.path.join(datadir,outRF))]
        outRF_new.write(pdbid + "," + lines[1].split(",")[1])
        outRF_new.close()
        os.system("mv " + os.path.join(datadir,"RF" + i + "_new.csv") + " " + os.path.join(datadir, "RF" + i + ".csv"))
         
        out.append(test_new)
        outRF = pd.read_csv(os.path.join(datadir,"RF" + i + ".csv"),dtype = {"pdb":str})
        out.append(outRF)


    os.chdir(datadir)
    get_output(out,outfile)
    os.system("rm " +  "RF*")
    os.chdir(olddir)



if __name__ == "__main__":

    main()



