"""
Vina 58 Features
"""
__author__ = "Jianing Lu"
__copyright__ = "Copyright 2018, NYU"
__license__ = ""

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import os
import fileinput
import sys
if sys.platform == "linux":
    import software_path_linux as path
elif sys.platform == "darwin":
    import software_path_mac as path




MGLPY = path.path_mgl_python()
MGLUTIL = path.path_mgl_script()
VINADIR = path.path_vina()
obable = path.path_obabel()
#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------
def runVina(fn,protpdbqt, ligpdbqt):
    """Run modified AutoDock Vina program with Vina score and 58 features """

    cmd = VINADIR + " --receptor " + protpdbqt + " --ligand " + ligpdbqt + "  --score_only --log score_v1.txt > out_vina.log"
    os.system(cmd)

    vinalist = [fn]
    for lines in fileinput.input("score_v1.txt"):
        if lines[0:4] in ["Affi", "Term"]:
            if sys.platform == "linux":
                line = lines.strip("\n").split()
                vinalist.append(line[1])
            elif sys.platform == "darwin":
                line = lines.strip("\n").split()
                vinalist.append(line[2])


    if len(vinalist) != 60:
        vinalist = [fn] + ['NA' for i in range(59)]

    return vinalist


def prepareProt(inprot, protpdbqt):
    """Prepare protein PDBQT file by MGLTools """

    cmd = MGLPY + " "  + MGLUTIL + "prepare_receptor4.py -r ../" + inprot + " -o " + protpdbqt + " -U 'nphs' > out1.tmp"
    os.system(cmd)


def prepareLig(inlig, ligpdbqt):
    """Prepare ligand PDBQT file by MGLTools """

    cmd = MGLPY + " "  + MGLUTIL + "prepare_ligand4.py -l ../" + inlig  + " -o " + ligpdbqt +  " -U 'nphs' > out2.tmp"
    os.system(cmd)

def featureVina(outfile, fn, inpro, inlig, datadir):
    """Get Vina score and Vina features """
    os.chdir(datadir)
    olddir = os.getcwd()
    os.system("mkdir Feature_Vina")
    os.chdir("Feature_Vina")
    protpdbqt = inpro.split(".")[0] + ".pdbqt"
    ligpdbqt = inlig.split(".")[0] + ".pdbqt"
    prepareProt(inpro,protpdbqt)
    prepareLig(inlig,ligpdbqt)
    vinalist = runVina(fn,protpdbqt,ligpdbqt)
    outfile.write( ",".join(vinalist) + "\n")
    os.chdir(olddir)
    os.system("rm -r Feature_Vina")

    return None

# for Vina flexible docking result, just do the score_only to get the Vina 58 features
def featureVina_flexible(outfile,fn, inpro, inlig, datadir):
    """Get Vina score and Vina features """
    os.chdir(datadir)
    olddir = os.getcwd()
    os.system("mkdir Feature_Vina")
    os.chdir("Feature_Vina")
    os.system("cp ../" + inlig + " .")
    protpdbqt = inpro.split(".")[0] + ".pdbqt"
    ligpdbqt = inlig
    prepareProt(inpro,protpdbqt)
    vinalist = runVina(fn,protpdbqt,ligpdbqt)
    outfile.write(",".join(vinalist) + "\n")
    os.chdir(olddir)
    os.system("rm -r Feature_Vina")

    return None


def main():
    args = sys.argv[1:]
    if args[-1] == "file":
        pdbfile = open('%s'%(sys.argv[1] + sys.argv[2]),'r')
        pdblist = []
        for i in pdbfile.readlines():
            pdblist.append(i[0:4])
    else:
        pdblist = []
        pdblist.append(sys.argv[2])
    datadir = sys.argv[1]
    outfile = open(datadir + "Vina58.csv","w")
    outfile.write('pdb,vina,' + ','.join(['vina' + str(n+1) for n in range(58)]) + "\n")
    outfile_RW = open(datadir + "Vina58_RW.csv","w")
    outfile_RW.write('pdb,vina,' + ','.join(['vina' + str(n+1) for n in range(58)]) + "\n")
    outfile_min = open(datadir + "Vina58_min.csv","w")
    outfile_min.write('pdb,vina,' + ','.join(['vina' + str(n+1) for n in range(58)]) + "\n")
    outfile_min_RW = open(datadir + "Vina58_min_RW.csv","w")
    outfile_min_RW.write('pdb,vina,' + ','.join(['vina' + str(n+1) for n in range(58)]) + "\n")
    for fn in pdblist:
        inpro = fn + "_protein.pdb"
        inpro_RW = fn + "_protein_RW.pdb"
        if fn + "_ligand.mol2" in os.listdir(datadir):
            inlig = fn + "_ligand.mol2"
        elif fn + "_ligand.sdf" in os.listdir(datadir):
            inlig =  fn + "_ligand.sdf"
            inlig_out =  fn + "_ligand.mol2"
            os.system(obable + " -isdf " + datadir + inlig + " -omol2 -O " + datadir + inlig_out)
            inlig = inlig_out
        elif fn + "_ligand_rigid.pdbqt" in os.listdir(datadir):
            inlig =  fn + "_ligand_rigid.pdbqt"
            inlig_out =  fn + "_ligand.mol2"
            os.system(obable + " -ipdbqt " + datadir + inlig + " -omol2 -O " + datadir + inlig_out)
            inlig = inlig_out
        elif fn + "_ligand_flexible.pdbqt" in os.listdir(datadir):
            inlig =  fn + "_ligand_flexible.pdbqt"
            inlig_out =  fn + "_ligand.mol2"
            os.system(obable + " -ipdbqt " + datadir + inlig + " -omol2 -O " + datadir + inlig_out)
            inlig = inlig_out

        else:
            print("wrong ligand input file format, it should be mol2, sdf, or pdbqt")
        featureVina(outfile,fn,inpro,inlig,datadir)
        featureVina(outfile_RW,fn,inpro_RW,inlig,datadir)
        inlig_min = fn + "_lig_min.pdb"
        featureVina(outfile_min,fn,inpro,inlig_min,datadir)
        inlig_min_RW = fn + "_lig_min_RW.pdb"
        featureVina(outfile_min_RW,fn,inpro_RW,inlig_min_RW, datadir)
        
    outfile.close()
    outfile_RW.close()
    outfile_min.close()
    outfile_min_RW.close()

    return None

if __name__ == "__main__":

    main()






