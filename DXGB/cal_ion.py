#-----------------------------------------------------------------------------
# Ion Feature
#-----------------------------------------------------------------------------
import os
import DXGB.get_pdbinfo as get_pdbinfo
from DXGB.get_pdbinfo import *
import numpy as np
import sys


def get_Ions(fn,lig,pro,infile):
    outfile = open(infile,"w")

    lig = get_pdbinfo.pdbinfo(fn,file = lig)
    lig_atoms = lig.getPolarAtoms()

    pro = get_pdbinfo.pdbinfo(fn,file = pro)
    pro_ions = pro.getIons()

    ion_coord = get_pdbinfo.pdbinfo(fn,lines = pro_ions).getCoords()
    lig_coord = get_pdbinfo.pdbinfo(fn,lines = lig_atoms).getCoords()
    
    ion_coord = np.expand_dims(ion_coord, 1)
    lig_coord = np.expand_dims(lig_coord,0)

    if ion_coord.shape[0] == 0:
        print("No Ion")
        outfile.close()
    else:
        distance = np.linalg.norm(ion_coord - lig_coord, axis = 2)
        distance_min = np.min(distance, axis = 1)

        ion_index = []
        for idx, i in enumerate(distance_min):
            if i < 3.5:
                ion_index.append(idx)

        for i in ion_index:
            ion_line = pro_ions[i]
            ion_distance = distance[i]
            ion_chain = get_pdbinfo.chid(ion_line)
            ion_name = get_pdbinfo.atmn(ion_line).strip()

            ion_chain = get_pdbinfo.chid(ion_line)
            if ion_chain != " ":
                ion_idx = str(int(get_pdbinfo.resi(ion_line))) + "." + ion_chain 
            else:
                ion_idx = str(int(get_pdbinfo.resi(ion_line)))
            
            for idx, d in enumerate(ion_distance):
                if d < 3.5:
                    lig_line = lig_atoms[idx]
                    lig_name = get_pdbinfo.atmn(lig_line).strip()
                    outline = fn + "," + ion_idx + "," + ion_name + "," + lig_name + "," + str(round(d,2)) + "\n"
                    outfile.write(outline)
        outfile.close()

def get_num(fn, infile):
    ion = []
    for line in open(infile):
        ion.append(line.split(",")[1] + line.split(",")[2])
    ionset = set(ion)
    return len(ionset)

def cal_Ni(outfile,fn, inprot, inlig, datadir):
    """
    Calculate Ion feature
    
    :param outfile: outfile object (can be wrote)
    :param fn: input index
    :param inprot: inpro
    :param inlig: inlig
    :param datadir: directory for input 

    return: the number of ions for input index will be wrote into the outfile

    """
    pro = os.path.join(datadir, inprot)
    lig = os.path.join(datadir,inlig)
    infile = os.path.join(datadir,"Ion_infor.dat")
    get_Ions(fn,lig,pro,infile)
    num = get_num(fn,infile)
    outfile.write(fn + "," + str(num) + "\n")

