#-----------------------------------------------------------------------------
# Vina Features
#-----------------------------------------------------------------------------

import os
import fileinput
import sys


def runVina(fn,protpdbqt, ligpdbqt):
    """Run modified AutoDock Vina program with Vina score and 58 features """

    cmd = "$VINADIR/vina --receptor " + protpdbqt + " --ligand " + ligpdbqt + "  --score_only --log score_v1.txt > out_vina.log"
    os.system(cmd)
    vinalist = [fn]
    for lines in fileinput.input("score_v1.txt"):
        if lines[0:4] in ["Affi", "Term"]:
            try:
                line = lines.strip("\n").split()
                vinalist.append(line[1])
            except:
                line = lines.strip("\n").split()
                vinalist.append(line[2])
    if len(vinalist) != 60:
        vinalist = [fn] + ['NA' for i in range(59)]
    return vinalist

def prepareProt(inprot, protpdbqt):
    """Prepare protein PDBQT file by MGLTools """
    cmd = "$MGLPY $MGLUTIL/prepare_receptor4.py -r ../" + inprot + " -o " + protpdbqt + " -U 'nphs' > out1.tmp"
    os.system(cmd)
    return None

def prepareLig(inlig, ligpdbqt):
    """Prepare ligand PDBQT file by MGLTools """
    cmd = "$MGLPY $MGLUTIL/prepare_ligand4.py -l ../" + inlig  + " -o " + ligpdbqt +  " -U 'nphs' > out2.tmp"
    os.system(cmd)
    return None

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







