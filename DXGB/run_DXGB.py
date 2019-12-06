#!/usr/bin/env python
import os
import sys
import click
import pandas as pd
from DXGB.convert_file import RF20_main as RF20_main
from DXGB.run_features import run_features as run_features
from DXGB.run_models import *

@click.command()
@click.option("--model", default="DXGB", show_default=True, help="model name")
@click.option("--modeldir",default="../Model", show_default=True, help="absolute model directory")
@click.option("--datadir", default="../Test", show_default=True, help="absolute data directory")
@click.option("--pdbid", default=None, show_default=True, help="pdbid, ligand input should be pdbid_ligand.mol2 or sdf,\nprotein input should be pdbid_protein_all.pdb")
@click.option("--outfile", default="score.csv",show_default=True, help="output filename")
@click.option("--runfeatures",is_flag=True, show_default=True, help="run features calculation")
@click.option("--water", default=None, show_default=True, help="water type")
@click.option("--opt", default=None, show_default=True, help="opt type")
@click.option("--rewrite", is_flag=True, help="rewrite protein part water generation, ligand optimization, ligand conformation generation or not")
@click.option("--average",is_flag=True, help="average for 10 models")
@click.option("--modelidx", default="1", show_default=True, help="model index")
@click.option("--featuretype", default="all", show_default=True, help="which feature will be calculated, options:all,Vina,SASA,BW,Ion,dE ")
@click.option("--runrf", is_flag=True, help="get deltaVinaRF20 scores")


def main(model, modeldir, datadir, pdbid, outfile, runfeatures, water, opt, rewrite, average, modelidx, featuretype, runrf):
    
    datadir = os.path.realpath(datadir)
    if pdbid:
        print("pdb index: " + pdbid  )
    print("file directory: " + datadir)
    
    if runfeatures:
        print("feature will be calculated:" + featuretype)
    
    print("output filename : " + outfile)
    olddir = os.getcwd()
    if runfeatures:
        if not water:
            water = "n"
        if not opt:
            opt = "n"
        run_features(datadir, pdbid, water_type = water, opt_type = opt, rewrite = rewrite, feature_type = featuretype)
        os.chdir(olddir)


    if opt == "rbwo":
        data_type = ["","_min","_min_RW","_min_BW"]
    elif opt == "rwo":
        data_type = ["_min_RW"]
    elif opt == "bwo":
        data_type == ["_min_BW"]
    elif opt == "pwo":
        data_type == ["_min_PW"]
    elif opt == "o":
        data_type = ["_min"]
    else:
        if water == "rbw":
            data_type = ["","_RW","_BW"]
        elif water == "rw":
            data_type = ["_RW"]
        elif water == "bw":
            data_type == ["_BW"]
        elif water == "pw":
            data_type == ["_PW"]
        else:
            data_type = [""]
    
    data_type_new = data_type
    out = []
    for idx, i in enumerate(data_type):
        modeldir = os.path.realpath(modeldir)
        inf = "Input" + data_type_new[idx] + ".csv"
        test_new = run_model(inf,datadir,i,model_dir = modeldir, model_name = model, average = average, model_index = modelidx)  
        out.append(test_new)
        if runrf:
            outRF = "RF" + data_type_new[idx] + ".csv"
            RF20_main(datadir,inf,outRF)
            ### correct wrong pdb name in scientific pattern ###
            outRF_new = open(os.path.join(datadir,"RF" + data_type_new[idx] + "_new.csv"),"w")
            outRF_new.write("pdb,RF20" + i + "\n")
       	    lines = [line for line in open(os.path.join(datadir,outRF))]
            if pdbid:
                outRF_new.write("".join([pdbid + "," + ",".join(line.split(",")[1:]) for line in lines[1:]]))
            else:
                pdbid = pd.read_csv(os.path.join(datadir,"Input" + data_type_new[idx] + ".csv"), dtype = {"pdb":str})["pdb"].tolist()
                outRF_new.write("".join([pdbid[idx] + "," + ",".join(line.split(",")[1:]) for idx, line in enumerate(lines[1:])]))
            outRF_new.close()
            os.system("mv " + os.path.join(datadir,"RF" + data_type_new[idx] + "_new.csv") + " " + os.path.join(datadir, "RF" + data_type_new[idx] + ".csv"))
         
            outRF = pd.read_csv(os.path.join(datadir,"RF" + data_type_new[idx] + ".csv"),dtype = {"pdb":str})
            out.append(outRF)

    os.chdir(datadir)
    get_output(out,outfile)

    if runrf:
        os.system("rm " +  "RF*")

    os.chdir(olddir)

if __name__ == "__main__":
    main()




