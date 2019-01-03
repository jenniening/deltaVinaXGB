import os
import sys
from combine_data import combine
from software_path import path_chimera
from software_path import path_python



chimera = path_chimera()
python = path_python()


def get_infor(datadir, pdbfile, features, type):
    if len(features) == 5:
        os.system(python + " cal_vina58.py " + datadir + " " + pdbfile + " " + type )
        print("Finish Vina58")
        os.system(python + " cal_sasa.py " + datadir + " " + pdbfile + " " + type)
        print("Finish SASA")
        os.system(python + " cal_dERMSD.py " + datadir + " " + pdbfile + " " + type)
        print("Finish ligand stability features")
        os.system(chimera + " --nogui --silent --script 'cal_bw.py " + datadir
                  + " " + pdbfile + " " + type + "'")
        print("Finish BW")
        os.system(chimera + " --nogui --silent --script 'cal_ion.py " + datadir
                  + " " + pdbfile + " " + type + "'")
        print("Finish Ion")
        combine(datadir)
        print("Combine all features data together")
    else:
        for f in features:
            if f == "vina58":
                os.system(python + " cal_vina58.py " + datadir + " " + pdbfile + " " + type )
                print("Finish Vina58")
            if f == "sasa":
                os.system(python + " cal_sasa.py " + datadir + " " + pdbfile + " " + type)
                print("Finish SASA")
            if f == "dE":
                os.system(python + " cal_RMSD.py " + datadir + " " + pdbfile + " " + type)
                print("Finish ligand stability features")
            if f == "bw":
                os.system(chimera + " --nogui --silent --script cal_bw.py " + datadir + " " + pdbfile + " " + type)
                print("Finish BW")
            if f == "ion":
                os.system(chimera + " --nogui --silent --script cal_ion.py " + datadir + " " + pdbfile + " " + type)
                print("Finish Ion")




def main():
    args = sys.argv[1:]
    if len(args) <= 2:
        print("usage: python run_all_features.py [--dir] dir [--PDB_id] PDB_id.txt [--f] features ")
        sys.exit(1)
    if not args:
        print("usage: python run_all_features.py [--dir] dir [--PDB_id] PDB_id.txt [--f] features ")
        sys.exit(1)
    elif sys.argv[1] == '--help':
        print("usage: python run_all_features.py [--dir] dir [--PDB_id] PDB_id.txt [--f] features ")
        sys.exit(1)
    if sys.argv[1] != "--dir":
        print("You should provide the directory of your data file by using --dir dir")
        sys.exit(1)
    else:
        datadir = sys.argv[2]


    if sys.argv[3] == "--PDB_id":

        if "." in sys.argv[4]:
            type = "file"
        else:
            type = "pdb"
        pdbfile = sys.argv[4]
    else:
        print("You should provide the pdb name of your data by using --PDB_id PDB_id.txt")
        sys.exit(1)
    if len(args) >= 6:
        if sys.argv[5] == "--f":
            if sys.argv[6] == "all":
                features = ["vina58","sasa","bw","ion","dE"]
            else:
                features = []
                features.append(sys.argv[6])
    else:
        features = ["vina58","sasa","bw","ion","dE"]
    print("pdb index file name: " + pdbfile)
    print("pdb index file type: " + type)
    print("file directory: " + datadir)
    print("features will be calcualted: " + ",".join(features))

    get_infor(datadir,pdbfile,features,type)

if __name__ == "__main__":

    main()



