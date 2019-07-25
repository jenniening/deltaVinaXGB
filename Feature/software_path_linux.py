# Change the softwares we need use (Linux)
import os
def path_python():
    python = "/beegfs/anaconda3_latest/bin/python"
    return python

def path_Rscript():
    R = "/share/apps/r/3.4.2/intel/bin/Rscript"
    return R

def path_obabel():
    obable = "/share/apps/openbabel/2.4.1/intel/bin/obabel"
    return obable

def path_mgl_python():
    MGLPY = "/share/apps/mgltools/1.5.6/bin/pythonsh"
    return MGLPY

def path_mgl_script():
    MGLUTIL = "/share/apps/mgltools/1.5.6/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/"
    return MGLUTIL

def path_vina():
    vina = os.path.realpath("../vina_package/bin/vina_linux")
    return vina

def path_msms():
    msmsdir = os.path.realpath("../msms_Linux_2.6.1/")
    return msmsdir

def path_RF20da():
    RF20da = os.path.realpath("RF20.rda")
    return RF20da
