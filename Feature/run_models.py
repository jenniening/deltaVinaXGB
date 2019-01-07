import numpy as np
import os
import pandas as pd
import pickle


def run_model(infile,datadir,rt, model_dir, model_name = "1", average = True, model_index = "1"):
    model_dir = os.path.realpath(model_dir)
    infile = os.path.join(datadir, infile)
    num = len([line for line in open(infile)]) - 1
    test_in = pd.read_csv(infile)

    # Feature names
    if model_name == "model_allfeatures":
        SASA_type = ["P","N","DA","D","A","AR","H","PL","HA","SA"]
        SASA = ["P2." + i for i in SASA_type] + ["P2dl." + i for i in SASA_type] + ["P2dp." + i for i in SASA_type]
        vina =  ["vina" + str(i) for i in range(1,59)]
        water = ["Nbw"] + ["Epw"] + ["Elw"]
        ion = ["Ni"]
        conformation =["dE_global","RMSD_global","number1"]
        featlist =  vina +  SASA  + ion + water +conformation
 
        model_dir = os.path.join(model_dir, model_name)

    elif model_name[0:14] == "model_fragment":
        SASA_type = ["P","N","DA","D","A","AR","H","PL","HA","SA"]
        SASA = ["P2." + i for i in SASA_type] + ["P2dp.SA","P2dl.SA"]
        vina =  ["vina_x","vina_y","vina1","vina3","vina53","vina55","vina54","vina56","vina4","vina52","vina58","vina48"]
        water = ["Nbw"] + ["Epw"] + ["Elw"]
        ion = ["Ni"]
        conformation =["dE_global","RMSD_global","number1"]
        featlist =  vina + SASA + ion + water +conformation + ["num_frag"]
        model_dir = os.path.join(model_dir,model_name)


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

    if rt == "":
        test_new.columns = ['pdb','vina','XGB']
    elif rt == "_min":
        test_new.columns = ['pdb','vina_min','XGB_min']
    elif rt == "_min_RW":
        test_new.columns = ['pdb','vina_min_RW','XGB_min_RW']
    elif rt == "_RW":
        test_new.columns = ['pdb','vina_RW','XGB_RW']

    return test_new

def get_output(out,outfile):
    n_out = len(out)
    df = out[0]
    for i in range(1,n_out):
        df = pd.merge(df,out[i],on = "pdb")
    df.to_csv(outfile, index = False)

    return None
