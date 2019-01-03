import os
import sys
import click
import pandas as pd

from software_path import path_python
from software_path import path_obabel
from convert_file import RF20_main
import run_features
from run_features import run_features
import run_models
from run_models import run_model
from run_models import get_output
import click




@click.command()
@click.option("--model", default = "model_allfeatures", show_default = True, help = "model name")
@click.option("--datadir", default = "/Users/jianinglu1/Documents/GitHub/deltaVinaXGB_develop/Test", show_default= True, help = "absolute data directory")
@click.option("--pdbid", default = "01", show_default = True, help = "pdbid, ligand input should be pdbid_ligand.mol2 or sdf,\nprotein input should be pdbid_protein_all.pdb")
@click.option("--outfile", default = "score.csv",show_default = True, help = "output filename")
@click.option("--runfeatures",is_flag = True, show_default = True, help = "run features calculation")
@click.option("--rw",is_flag = True, help = "receptor water flag")
@click.option("--water", is_flag = True, help = "water flag")
@click.option("--decoy", is_flag = True, help = "decoy flag")
@click.option("--average",is_flag = True, help = "average for 10 models")
@click.option("--modelidx", default = "1", show_default = True, help = "model index")

def main(model, datadir, pdbid, outfile, runfeatures, rw, water, decoy, average, modelidx):
    print("pdb index: " + pdbid  )
    print("file directory: " + datadir)
    print("output filename : " + outfile)
    olddir = os.getcwd()
    if runfeatures:
        run_features(datadir, pdbid, rw = rw, water = water, decoy = decoy)
        os.chdir(olddir)

    if water:
        data_type = ["","_min","_min_RW"]
    else:
        data_type = ["","_min"]

    if decoy:
        pass
    
    out = []
    for i in data_type:
        inf = "Input" + i + ".csv"
        test_new = run_model(inf,datadir,i,model_name = model, average = average, model_index = modelidx)
        outRF = "RF" + i + ".csv"
        RF20_main(datadir,inf,outRF)
        outRF = pd.read_csv(os.path.join(datadir,outRF))
        outRF.columns = ["pdb","RF20" + i]
        out.append(test_new)
        out.append(outRF)


    os.chdir(datadir)
    get_output(out,outfile)
    os.system("rm " +  "RF*")
    os.chdir(olddir)



if __name__ == "__main__":

    main()



