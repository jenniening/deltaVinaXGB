import numpy as np
import os
import pandas as pd
import pickle

def run_model(infile, datadir, rt, model_dir, model_name="DXGB", average=True, model_index="1"):
    """
    run predictions by our model
    
    :param infile: input file
    :param datadir: directory for input file
    :param rt: input structure type, can be "", "_min", "_min_RW", "_RW", "_min_BW", "_min_PW"
    :param model_dir: directory for models
    :param model_name: model name, defaults to "DXGB"
    :param average: whether to use ensemble model prediction, defaults to True
    :param model_index: if not average, which model should be used, defaults to "1"
    :return: test_new, the dataframe includes pdb, vina, and XGB scores.

    """

    model_dir = os.path.realpath(model_dir)
    infile = os.path.join(datadir, infile)
    num = len([line for line in open(infile)]) - 1
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

    test_new = pd.concat([test[["pdb","vina"]],test_pred],axis = 1) 
    test_new.columns = ["pdb", "vina" + rt, "XGB" + rt]

    return test_new

def get_output(out,outfile):
    """
    Merge output of XGB and RF20
    
    :param out: outputs are needed to merge
    :param outfile: filename for the merged results

    """
    n_out = len(out)
    df = out[0]
    for i in range(1,n_out):
        df = pd.merge(df,out[i],on = "pdb")
    df.to_csv(outfile, index = False)

