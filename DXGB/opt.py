#-----------------------------------------------------------------------------
# Optimized Structure
#-----------------------------------------------------------------------------
import os
import fileinput
import sys
import DXGB.get_pdbinfo as get_pdbinfo
from DXGB.utils import get_tool

vina = get_tool("vina")
prepare_ligand4 = get_tool("prepare_ligand4.py")
prepare_receptor4 = get_tool("prepare_receptor4.py")
obabel = get_tool("obabel")

#### Do minimization, structure file should have hydrogen atom!!!####
## generate the refined data list
## generate the box for mol2 file
def get_box(fn, inlig):
    if inlig.split(".")[-1] == "mol2":
        inputfile = open("../" + inlig)
        x = []
        y = []
        z = []
        flag = False
        for line in inputfile:
            if line[0:13] == "@<TRIPOS>ATOM":
                flag = True
                continue
            elif line[0:13] == "@<TRIPOS>BOND":
                flag = False
            if flag:
                x.append(float(line[16:26].split()[0]))
                y.append(float(line[27:37].split()[0]))
                if line[37] == "0":
                    z.append(float(line[38:48].split()[0]))
                else:
                    z.append(float(line[37:48].split()[0]))
    elif inlig.split(".")[-1] == "pdb":
        lines = get_pdbinfo.pdbinfo(name = fn, file = "../" + inlig).getAtoms()
        x,y,z = [],[],[]
        for line in lines:
            coords = get_pdbinfo.pdbinfo(name = fn, lines = [line]).getCoords()
            x.append(float(coords[0][0]))
            y.append(float(coords[0][1]))
            z.append(float(coords[0][2]))

    x_center =(max(x)+min(x))/2
    y_center =(max(y)+min(y))/2
    z_center = (max(z)+min(z))/2
    size_x = max(x) - min(x) + 10
    size_y = max(y) - min(y) + 10
    size_z = max(z) - min(z) + 10
    new_file = open("box.txt","w")
    new_file.write("center_x = " + str(x_center)+ "\n")
    new_file.write("center_y = " + str(y_center)+ "\n")
    new_file.write("center_z = " + str(z_center) + "\n")

    new_file.write("size_x = " + str(size_x)+ "\n")
    new_file.write("size_y = " + str(size_y)+ "\n")
    new_file.write("size_z = " + str(size_z) + "\n")
    new_file.close()

def genpdbqt(fn, ligpdb, propdb):
    propdbqt = fn + "_rec.pdbqt"
    ligpdbqt = fn + "_lig.pdbqt"
    cmd1 =  prepare_receptor4 + " -r ../"  + propdb + " -o " + propdbqt + " -U '_' > out1.tmp"
    cmd2 =  prepare_ligand4 + " -l ../" + ligpdb  + " -o " + ligpdbqt +  " -U '_' -Z > out2.tmp"
    os.system(cmd1)
    os.system(cmd2)


def runmin(fn):
    cmd = vina + " --receptor " + fn + "_rec.pdbqt --ligand " + fn + "_lig.pdbqt --config box.txt --local_only --out " + fn + "_lig_min.pdbqt >out_min.txt"
    os.system(cmd)


def chanPdb(fn):
    cmd = obabel + " -ipdbqt " + fn + "_lig_min.pdbqt  -opdb -O " + fn + "_lig_min.pdb"
    os.system(cmd)


def get_Co(datadir,fn, inlig, st):
    """
    Get Vina optimized structure 

    :param datadir:datadir for input 
    :param fn: input index
    :param inlig: inlig 
    :param st: water type

    """
    os.chdir(datadir)
    olddir = os.getcwd()
    os.system("mkdir vinamin_rigid")
    os.chdir("vinamin_rigid")
    if st != "" and "_" not in st:
        st = "_" + st
    
    inpro = fn + "_protein" + st + ".pdb"
    outlig = fn + "_lig_min" + st + ".pdb"
    
    genpdbqt(fn,inlig,inpro)
    get_box(fn,inlig)
    runmin(fn)
    chanPdb(fn)
    os.system("cp " + fn + "_lig_min.pdb ../" + outlig)
    os.chdir(olddir)
    os.system("rm -r vinamin_rigid")






