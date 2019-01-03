import os
import sys
try:
    sys.path.remove('/share/apps/chimera/1.11.2/lib/python2.7/site-packages')
except ValueError:
    pass
from combine_data_addCw import combine
import pandas as pd
import numpy as np
import pickle
from software_path import path_chimera
from software_path import path_python
from Convert_File import RF20_main
sys.path.append(os.path.dirname(os.path.realpath(__file__)))


chimera = path_chimera()
python = path_python()


def get_infor(datadir, pdbfile, features, file_type):
    # no remove
    os.system(chimera + " --nogui --silent --script 'get_Crw_noremove.py " + datadir + " " + pdbfile + " " + file_type + "'")
    #os.system(chimera + " --nogui --silent --script 'get_Crw.py " + datadir + " " + pdbfile + " " + file_type + "'")
    print("Finish generate RW")
    os.system(chimera + " --nogui --silent --script 'get_Co.py " + datadir + " " + pdbfile + " " + file_type + "'" )
    print("Finish generate Co, Corw")
    if len(features) == 5:
        os.system(python + " cal_vina58_addCw.py " + datadir + " " + pdbfile + " " + file_type )
        print("Finish Vina")
        os.system(python + " cal_sasa_addCw.py " + datadir + " " + pdbfile + " " + file_type)
        print("Finish SASA")
        os.system(python + " cal_dERMSD_correct.py " + datadir + " " + pdbfile + " " + file_type)
        print("Finish ligand stability features")
        # obabel water
        os.system(chimera + " --nogui --silent --script 'cal_bw_addCw.py " + datadir + " " + pdbfile + " " + file_type + "'")

        #os.system(chimera + " --nogui --silent --script 'cal_bw.py " + datadir + " " + pdbfile + " " + file_type + "'")
        print("Finish BW")
        os.system(chimera + " --nogui --silent --script 'cal_ion_addCw.py " + datadir + " " + pdbfile + " " + file_type + "'")
        print("Finish Ion")
        combine(datadir)
        print("Combine all features data together")
    else:
        for f in features:
            if f == "vina":
                os.system(python + " cal_vina58_addCw.py " + datadir + " " + pdbfile + " " + file_type )
                print("Finish Vina")
            if f == "sasa":
                os.system(python + " cal_sasa_addCw.py " + datadir + " " + pdbfile + " " + file_type)
                print("Finish SASA")
            #if f == "dE":
                #os.system(python + " cal_dERMSD.py " + datadir + " " + pdbfile + " " + file_type)
                #print("Finish ligand stability features")
            #### add by Yuwei 
            if f == "dE":
                os.system(python + " cal_dERMSD_correct.py " + datadir + " " + pdbfile + " " + file_type)
                print("Finish ligand stability features")
            ##### add end 

            if f == "bw":
                os.system(chimera + " --nogui --silent --script 'cal_bw_addCw.py " + datadir + " " + pdbfile + " " + file_type + "'")
                print("Finish BW")
            if f == "ion":
                os.system(chimera + " --nogui --silent --script 'cal_ion_addCw.py " + datadir + " " + pdbfile + " " + file_type + "'")
                print("Finish Ion")
       
        #######################
        # add by Yuwei Y 
        combine(datadir)
        print('Combine all features data together')

        # add end 
        #######################3


def run_model(infile,datadir,rt):
    number = sum(1 for line in open(datadir + infile)) -1
    test_in = pd.read_csv(datadir + infile)

    # Feature names
    SASA = ["P2.P","P2.N","P2.DA","P2.D","P2.A","P2.AR","P2.H","P2.PL","P2.HA","P2.SA"]
    SASA_dl = ["P2dl.P","P2dl.N","P2dl.DA","P2dl.D","P2dl.A","P2dl.AR","P2dl.H","P2dl.PL","P2dl.HA","P2dl.SA"]
    SASA_dp = ["P2dp.P","P2dp.N","P2dp.DA","P2dp.D","P2dp.A","P2dp.AR","P2dp.H","P2dp.PL","P2dp.HA","P2dp.SA"]
    vina10 =  ["vina1","vina3","vina53","vina55","vina54","vina56","vina4","vina52","vina58","vina48"]
    vina = ["vina" + str(i) for i in range(1,59)]
    water = ["Nbw"] + ["Epw"] + ["Elw"]
    ion = ["Ni"]
    conformation =["dE"] + ["RMSD"] + ["number1"]
    featlist =  vina +  SASA + SASA_dl + SASA_dp + ion + water +conformation 


    test = test_in.copy()
    test_features = test[featlist].values

    y_test = pd.DataFrame(np.zeros((number,10)),columns = list(str(i) + "_t" for i in range(1,11)))
    for i in range(1,11):
        mod = pickle.load(open("../Model/model_latest/pima.pickle_" + str(i) +  ".dat", "rb"))
        test_preds = mod.predict(test_features)
        y_test[str(i) + "_t"]= test_preds
    test_pred = y_test.sum(axis = 1)/10 + test.vina
    test_new = pd.concat([test[["pdb","vina"]],test_pred],axis = 1) 
    if rt == "C":
        test_new.columns = ['pdb','vina','XGB']
    elif rt == "Co":
        test_new.columns = ['pdb','vina_min','XGB_min']
    elif rt == "Corw":
        test_new.columns = ['pdb','vina_min_RW','XGB_min_RW']
    elif rt == "Crw":
        test_new.columns = ['pdb','vina_RW','XGB_RW']

    #test_new.to_csv(datadir + outfile,index = False)

    return test_new

def get_output(out, datadir,outfile):
    n_out = len(out)
    df = out[0]
    for i in range(1,n_out):
        df = pd.merge(df,out[i],on = "pdb")
    df.to_csv(datadir + outfile, index = False)

    return None


def main():
    args = sys.argv[1:]
    if len(args) <= 2:
        print("usage: python run_DXGB_directly.py [--dir] dir [--PDB_id] PDB_id.txt [--o] outfile")
        sys.exit(1)
    if not args:
        print("usage: python run_DXGB_directly.py [--dir] dir [--PDB_id] PDB_id.txt [--o] outfile")
        sys.exit(1)
    elif sys.argv[1] == '--help':
        print("usage: python run_DXGB_directly.py [--dir] dir [--PDB_id] PDB_id.txt [--o] outfile")
        sys.exit(1)
    if sys.argv[1] != "--dir":
        print("You should provide the directory of your data file by using --dir dir")
        sys.exit(1)
    else:
        if sys.argv[2].split("/")[0] == "..":
            current_dir = os.getcwd().split("/")[0:-1]
            datadir ="/".join(current_dir) + "/" + "/".join(sys.argv[2].split("/")[1:])
        else:
            datadir = sys.argv[2]
        if datadir[-1] != "/":
            datadir = datadir + "/"

    if sys.argv[3] == "--PDB_id":

        if "." in sys.argv[4]:
            file_type = "file"
        else:
            file_type = "pdb"
        pdbfile = sys.argv[4]
    else:
        print("You should provide the pdb name of your data by using --PDB_id PDB_id.txt")
        sys.exit(1)
    if sys.argv[5] == "--o":
        outfile = sys.argv[6]
    else:
        print("You should provide inputfile name by using --o outfile")
        sys.exit(1)

    features = ["vina","sasa","bw","ion","dE"]

    print("pdb index file name: " + pdbfile)
    print("pdb index file type: " + file_type)
    print("file directory: " + datadir)
#    print("features will be calcualted: " + ",".join(features))
    print("outfile name : " + outfile)
    get_infor(datadir,pdbfile,features,file_type)
    inf1 = "Input.csv"
    test_new1 = run_model(inf1,datadir,rt = "C")
    outRF1 = "RF1.csv"
    RF20_main(datadir,inf1,outRF1)
    inf2 = "Input_min.csv"
    test_new2 = run_model(inf2,datadir,rt = "Co")
    outRF2 = "RF2.csv"
    RF20_main(datadir,inf2,outRF2)
    inf3 = "Input_min_RW.csv"
    test_new3 = run_model(inf3,datadir,rt = "Corw")
    outRF3 = "RF3.csv"
    RF20_main(datadir,inf3,outRF3)
    inf4 = "Input_RW.csv"
    test_new4 = run_model(inf4, datadir,rt = "Crw")
    outRF4 = "RF4.csv"
    RF20_main(datadir, inf4, outRF4)
    outRF1 = pd.read_csv(datadir + outRF1)
    outRF1.columns = ['pdb','RF20']
    outRF2 = pd.read_csv(datadir + outRF2)
    outRF2.columns = ['pdb','RF20_min']
    outRF3 = pd.read_csv(datadir + outRF3)
    outRF3.columns = ['pdb','RF20_min_RW']
    outRF4 = pd.read_csv(datadir + outRF4)
    outRF4.columns = ['pdb','RF20_RW']
    get_output([test_new1,outRF1,test_new2,outRF2, test_new3,outRF3, test_new4,outRF4],datadir, outfile)
    os.system("rm " + datadir + "RF*")


if __name__ == "__main__":

    main()



