# Change the softwares we need use (Mac)
import os
def path_python():
    python = "/Users/jianinglu1/anaconda3/envs/DXG/bin/python"
    return python

def path_Rscript():
    R = "/usr/local/bin/Rscript"
    return R

def path_obabel():
    obable = "/Users/jianinglu1/anaconda3/envs/DXG/bin/obabel"
    return obable

def path_mgl_python():
    MGLPY = "/Users/jianinglu1/Documents/mgltools_i86Darwin9_1.5.6/bin/pythonsh"
    return MGLPY

def path_mgl_script():
    MGLUTIL = "/Users/jianinglu1/Documents/mgltools_i86Darwin9_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/"
    return MGLUTIL

def path_vina():
    vina = os.path.realpath("../vina_package/bin/vina_mac")
    return vina

def path_msms():
    msmsdir = os.path.realpath("../msms_MacOSX_2.6.1/")
    return msmsdir

def path_RF20da():
    RF20da = os.path.realpath("RF20_rm2016.rda")
    return RF20da
