#!/usr/bin/env python
import os
import sys
import pandas as pd
from DXGB.convert_file import RF20_main as RF20_main
from DXGB.run_features import run_features as run_features
from DXGB.run_models import run_model
from DXGB.run_models import get_output



def get_DXGB(model, modeldir, datadir, pdbid, outfile, runfeatures, water, opt, rewrite, average=True, modelidx="1", featuretype="all", runrf=False, runscore=True):
    """
    Get deltaVinaXGB all features and score
    
    Parameters
    ----------
    model : str
        model type, "DXGB" is the one for our previous trained deltaVinaXGB.
    modeldir : str
        directory for models, the previous trained models should be in the modeldir/model.
    datadir : str
        directory for data(input structures or features), and all output files.
    pdbid : str
        unique index for input data, can be pdb id, or any other customrized index.
    outfile : str
        output score file name (format is csv)
    runfeatures : bool
        whether to calculate features.
    water : str
        water type, can be:
        "rbw" --> all types of water, 
        "rw" --> only receptor water, 
        "bw" --> only bridging water,
        False --> no water.
    opt : str
        optimization type, can be:
        "rbwo" --> all types of optimization, 
        "rwo" --> only optimizaion with receptor water,
        "bwo" --> only optimization with bridging water, 
        "o" --> optimization with no water,
        False --> no optimization.
    rewrite : bool
        whether to rewrite protein with water stucture, optimized ligand structure, ligand conformation generation results.
    average : bool
        whether to use ensemble model results, defaults to True.
    modelidx : str
        requires when average=False, used to define which model will be used to predict score, defaults to "1".
    featuretype : str
        requires when runfeatures=True, used to define which features will be calculated, defaults to "all", can be:
        "all" --> calculate all features,
        "Vina" --> calculate Vina58, 
        "SASA" --> calculate SASA, 
        "BW" --> calculate bridging water features, 
        "Ion" --> calculate ion features, 
        "dE" --> calculate ligand stability features.
    runrf : bool
        whether to calculate deltaVinaRF score, defaults to False.
    runscore:bool
        whether to calculate score, defaults to True. If False and runfeatures=True, only performes feature calculation.
        
    Return
    ----------
    output score will be saved in outfile, together with all other calculated features will be saved in datadir

    """
    
    datadir = os.path.realpath(datadir)
    if pdbid:
        print("pdb index: " + pdbid  )
    print("file directory: " + datadir)
    
    if runfeatures:
        print("feature will be calculated:" + featuretype)
    
    olddir = os.getcwd()
    if runfeatures:
        if not water:
            water = "n"
        if not opt:
            opt = "n"
        run_features(datadir, pdbid, water_type = water, opt_type = opt, rewrite = rewrite, feature_type = featuretype)
        os.chdir(olddir)

    if runscore:
        print("output score file: " + outfile)
        modeldir = os.path.realpath(modeldir)

        if opt == "rbwo":
            data_type = ["","_min","_min_RW","_min_BW"]
        elif opt == "rwo":
            data_type = ["_min_RW"]
        elif opt == "bwo":
            data_type == ["_min_BW"]
        elif opt == "o":
            data_type = ["_min"]
        else:
            if water == "rbw":
                data_type = ["","_RW","_BW"]
            elif water == "rw":
                data_type = ["_RW"]
            elif water == "bw":
                data_type == ["_BW"]
            else:
                data_type = [""]
    
        data_type_new = data_type
        out = []
        for idx, i in enumerate(data_type):
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





