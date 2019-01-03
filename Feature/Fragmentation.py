import os

fragment = "/scratch/jl7003/BACE/script/CyclicPoseFragmentation_noH.py"

datadir = "/scratch/jl7003/BACE/prepare_final/"
pdblist = ['{:0>3}'.format(i) for i in range(1,155)]
olddir = os.getcwd()
for fn in pdblist:
    os.chdir(datadir)
    os.system("mkdir pose_final")
    os.chdir("pose_final")
    os.system("mkdir " + fn )
    os.chdir(fn)
    os.system("cp ../../"  + fn +  "/" + fn + "_ligand.sdf .")
    in_lig = fn + "_ligand.sdf"
    os.system("python " + fragment + " --input " + in_lig)
    os.chdir(olddir)
