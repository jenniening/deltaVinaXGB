import numpy as np
import os
import pandas as pd
import pickle

def run_model(infile, datadir,rt, model_dir, model_name="DXGB", average=True, model_index="1", decoy=False):
    model_dir = os.path.realpath(model_dir)
    infile = os.path.join(datadir, infile)
    num = len([line for line in open(infile)]) - 1

    if decoy:
        test_in = pd.read_csv(infile,dtype = {"pdb":str,"idx":str})
    else:
        test_in = pd.read_csv(infile,dtype = {"pdb":str})

    # Feature names
    model_dir = os.path.join(model_dir, model_name)
    featlist = [line.rstrip() for line in open(os.path.join(model_dir, "featlist.csv"))]

    test = test_in.copy()
    test_features = test[featlist].values

    if average:
        y_test = pd.DataFrame(np.zeros((num,10)),columns = list(str(i) + "_t" for i in range(1,11)))
        for i in range(1,11):
            print(i)
            mod = pickle.load(open(model_dir + "/pima.pickle_" + str(i) +  ".dat", "rb"))
            test_preds = mod.predict(test_features)
            y_test[str(i) + "_t"]= test_preds
        test_pred = y_test.sum(axis = 1)/10 + test.vina
    else:
        y_test = pd.DataFrame(np.zeros((num,1)),columns = list(str(i) + "_t" for i in range(1,2)))
        mod = pickle.load(open(model_dir + "/pima.pickle_" + str(model_index) +  ".dat", "rb"))
        test_preds = mod.predict(test_features)
        y_test["1_t"]= test_preds
        test_pred = y_test.sum(axis = 1)/1 + test.vina
    if decoy:
        test_new = pd.concat([test[["pdb","idx","vina"]],test_pred],axis = 1)
    else:
        test_new = pd.concat([test[["pdb","vina"]],test_pred],axis = 1) 
    
    if decoy:
        if rt == "":
            test_new.columns = ['pdb','idx','vina','XGB']
        elif rt == "_min":
            test_new.columns = ['pdb','idx','vina_min','XGB_min']
        elif rt == "_min_RW":
            test_new.columns = ['pdb','idx','vina_min_RW','XGB_min_RW']
        elif rt == "_RW":
            test_new.columns = ['pdb','idx','vina_RW','XGB_RW']
    else:

        if rt == "":
            test_new.columns = ['pdb','vina','XGB']
        elif rt == "_min":
            test_new.columns = ['pdb','vina_min','XGB_min']
        elif rt == "_min_RW":
            test_new.columns = ['pdb','vina_min_RW','XGB_min_RW']
        elif rt == "_RW":
            test_new.columns = ['pdb','vina_RW','XGB_RW']
        elif rt == "_min_BW":
            test_new.columns = ['pdb','vina_min_BW','XGB_min_BW']
        elif rt == "_min_PW":
            test_new.columns = ['pdb','vina_min_PW','XGB_min_PW']

    return test_new

def get_output(out,outfile,decoy):
    n_out = len(out)
    df = out[0]
    for i in range(1,n_out):
        if decoy:
            df = pd.merge(df,out[i],on = ["pdb","idx"])
        else:
            df = pd.merge(df,out[i],on = "pdb")
    df.to_csv(outfile, index = False)

    return None
