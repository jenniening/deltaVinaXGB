import os

# Change the softwaters we need use here
def path_chimera():
    chimera ="/share/apps/chimera/1.11.2/bin/chimera"
    return chimera

def path_python():
    python = "~/anaconda3/bin/python"
    return python

def path_mgl_python():
    MGLPY = "/share/apps/mgltools/1.5.6/bin/pythonsh"
#    "/Users/jianinglu1/Documents/mgltools_i86Darwin9_1.5.6/bin/pythonsh"
    return MGLPY

def path_mgl_script():
    MGLUTIL = "/share/apps/mgltools/1.5.6/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/"
#    "/Users/jianinglu1/Documents/mgltools_i86Darwin9_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/"
    return MGLUTIL

def path_vina():
    vina = "/scratch/jl7003/deltaXGB/vina_package/bin/vina_linux"
    return vina

def path_obabel():
    obable = "/share/apps/openbabel/2.4.1/intel/bin/obabel"
    return obable

def path_msms():
    msmsdir = "/home/jl7003/msms/"
    return msmsdir

def path_R():
    R = "/share/apps/r/3.4.2/intel/bin/R"
    return R
