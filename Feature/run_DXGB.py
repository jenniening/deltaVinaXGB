import os
import sys
import click
import pandas as pd
from Feature.convert_file import RF20_main as RF20_main
from Feature.run_features import run_features as run_features
from Feature.run_models import *

@click.command()
@click.option("--model", default="DXGB", show_default=True, help="model name")
@click.option("--modeldir",default="../Model", show_default=True, help="absolute model directory")
@click.option("--datadir", default="../Test", show_default=True, help="absolute data directory")
@click.option("--pdbid", default=None, show_default=True, help="pdbid, ligand input should be pdbid_ligand.mol2 or sdf,\nprotein input should be pdbid_protein_all.pdb")
@click.option("--outfile", default="score.csv",show_default=True, help="output filename")
@click.option("--decoy", is_flag=True, help="decoy flag")
@click.option("--decoydatadir", default=None, show_default=True, help="decoy datadir, if decoy == True, please provide decoydatadir, and datadir is the directory for reference file")
@click.option("--prodatadir",default=None, show_default=True, help="protein directory, needed when protein and ligand are not in same directory")
@click.option("--proid", default=None, show_default=True, help="protein id, needed when protein and ligand have different ids")
@click.option("--runfeatures",is_flag=True, show_default=True, help="run features calculation")
@click.option("--water", default=None, show_default=True, help="water type")
@click.option("--opt", default=None, show_default=True, help="opt type")
@click.option("--decoytype",default="docking",show_default=True, help="decoy type")
@click.option("--rewrite", is_flag=True, help="rewrite protein part water generation, ligand optimization, ligand conformation generation or not")
@click.option("--average",is_flag=True, help="average for 10 models")
@click.option("--modelidx", default="1", show_default=True, help="model index")
@click.option("--ligname", is_flag=True, help="whether use pdbid to get decoys with same name (CASF-2013/2016 screening)")
@click.option("--featuretype", default="all", show_default=True, help="which feature will be calculated, options:all,Vina,SASA,BW,Ion,dE ")
@click.option("--runrf", is_flag=True, help="get deltaVinaRF20 scores")


def main(model, modeldir, datadir, decoydatadir, prodatadir, pdbid, proid, outfile, runfeatures, water, opt, decoy, decoytype, rewrite, average, modelidx,ligname, featuretype, runrf):
    """
    :param datadir: directory for input structures, files and output scores
    :param pdbid: used to find the corrsponding structure files 
    :param outfile: outfile name, defaults to score.csv
    :param runfeatures: defaults to all, can be "all", "Vina", "SASA", "BW", "Ion", "dE"
    :param water: defaults to None, can be "rbw", "rw", "bw". To get all types of water, use "rbw"
    :param opt: defaults to None, can be "rbwo", "rwo", "bwo", "o". To get all types of optimizations, use "rbwo"
    """
    
    datadir = os.path.realpath(datadir)
    if pdbid:
        print("pdb index: " + pdbid  )
    print("file directory: " + datadir)
    if decoy:
        decoydatadir = os.path.realpath(decoydatadir)
        print("pdb index: " + pdbid  )
        print("ref directory: " + datadir)
        print("decoy directory: " + decoydatadir)
        print("decoy type: " + decoytype)
    if prodatadir:
        print("protein directory is not same as ref directory: " + prodatadir)
    if proid:
        print("protein id is not same as pdb index:" + proid)
    
    if runfeatures:
        print("feature will be calculated:" + featuretype)
    
    print("output filename : " + outfile)
    olddir = os.getcwd()
    if runfeatures:
        if not water:
            water = "n"
        if not opt:
            opt = "n"
        run_features(datadir, prodatadir, decoydatadir, pdbid, proid, water_type = water, opt_type = opt, decoy_type = decoytype, rewrite = rewrite, decoy = decoy, ligname = ligname, feature_type = featuretype)
        os.chdir(olddir)

    if decoy:
        if opt == "rbwo":
            data_type = ["","_min","_min_RW","_min_BW"]
        elif opt == "rwo":
            data_type = ["_min_RW"]
        elif opt == "bwo":
            data_type = ["_min_BW"]
        elif opt == "pwo":
            data_type = ["_min_PW"]
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
    else:

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
    
    if decoy:
        datadir = decoydatadir
        data_type_new = ["_decoys" + i for i in data_type]
    else:
        data_type_new = data_type
    out = []
    for idx, i in enumerate(data_type):
        inf = "Input" + data_type_new[idx] + ".csv"
        test_new = run_model(inf,datadir,i,model_dir = modeldir, model_name = model, average = average, model_index = modelidx, decoy = decoy)  
        out.append(test_new)
        if runrf:
            outRF = "RF" + data_type_new[idx] + ".csv"
            RF20_main(datadir,inf,outRF, decoy)
            ### correct wrong pdb name in scientific pattern ###
            outRF_new = open(os.path.join(datadir,"RF" + data_type_new[idx] + "_new.csv"),"w")
            if decoy:
                outRF_new.write("pdb,idx,RF20" + i + "\n")
            else:
                outRF_new.write("pdb,RF20" + i + "\n")
       	    lines = [line for line in open(os.path.join(datadir,outRF))]
            if pdbid:
                outRF_new.write("".join([pdbid + "," + ",".join(line.split(",")[1:]) for line in lines[1:]]))
            else:
                pdbid = pd.read_csv(os.path.join(datadir,"Input" + data_type_new[idx] + ".csv"), dtype = {"pdb":str})["pdb"].tolist()
                outRF_new.write("".join([pdbid[idx] + "," + ",".join(line.split(",")[1:]) for idx, line in enumerate(lines[1:])]))
            outRF_new.close()
            os.system("mv " + os.path.join(datadir,"RF" + data_type_new[idx] + "_new.csv") + " " + os.path.join(datadir, "RF" + data_type_new[idx] + ".csv"))
         
            if decoy:
                outRF = pd.read_csv(os.path.join(datadir,"RF" + data_type_new[idx] + ".csv"),dtype = {"pdb":str,"idx":str})
            else:
                outRF = pd.read_csv(os.path.join(datadir,"RF" + data_type_new[idx] + ".csv"),dtype = {"pdb":str})
            out.append(outRF)

    os.chdir(datadir)
    get_output(out,outfile,decoy)
    if runrf:
        os.system("rm " +  "RF*")
    os.chdir(olddir)

    return None

if __name__ == "__main__":
    main()



