"""
Pharamcophore based SASA for complex
"""

__author__ = "Jianing Lu"
__copyright__ = "Copyright 2018, NYU"
__license__ = ""

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import os, sys, pybel
try: 
    sys.path.remove('/share/apps/chimera/1.11.2/lib/python2.7/site-packages')
except ValueError:
    pass
#print(sys.path)
import openbabel as ob
import numpy as np
import pandas as pd

from pharma import pharma
from software_path import path_msms

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------

msmsdir = path_msms()

def runMSMS(inprot, inlig, MSMSDIR = '.'):
    olddir = os.getcwd()
    os.chdir(MSMSDIR)
    os.system('mkdir tmp_dl')
    ppdb = 'tmp_dl/p.pdb'
    __, intype = os.path.splitext(inprot)
    if intype[1:].lower() != 'pdb':
        prot = pybel.readfile(intype[1:], inprot).__next__()
        output = pybel.Outputfile("pdb", ppdb, overwrite=True)
        output.write(prot)
        output.close()
    else:
        os.system("""sed 's/HETATM/ATOM\ \ /g' """ + inprot + " > " + ppdb)
    lpdb = 'tmp_dl/l.pdb'
    __, intype = os.path.splitext(inlig)
    if intype[1:].lower() != 'pdb':
        lig = pybel.readfile(intype[1:], inlig).__next__()
        output = pybel.Outputfile("pdb", lpdb, overwrite=True)
        output.write(lig)
        output.close()
    else:
        os.system("""sed 's/HETATM/ATOM\ \ /g' """ + inlig + " > " + lpdb)

    os.chdir('tmp_dl')
    os.system("cp " + msmsdir + "atmtypenumbers .")
    ppdb2 = 'p_sa.pdb'
    lpdb2 = 'l_sa.pdb'

    pidx, ppharm = pharma('p.pdb').assign(write=True, outfn = ppdb2)
    lidx, lpharm = pharma('l.pdb').assign(write=True, outfn = lpdb2)
    elementint = [6, 7, 8, 9, 15, 16, 17, 35, 53]
    psub = [idx for idx in pidx if ppharm[idx][0] in elementint]
    lsub = [idx for idx in lidx if lpharm[idx][0] in elementint]
    comp = []
    for idx in lsub:
        comp.append(lpharm[idx][0:2])

    df1 = {}
    df1['atm'] = np.array(comp)[:,0]
    df1['pharma'] = np.array(comp)[:,1]
    df1 = pd.DataFrame(df1)
    # pdb to xyzr convert
    os.system(msmsdir + "pdb_to_xyzr " + ppdb2 + " > p_sa.xyzr")
    os.system(msmsdir + "pdb_to_xyzr " + lpdb2 + " > l_sa.xyzr")
    os.system("cat p_sa.xyzr l_sa.xyzr > pl_sa.xyzr")

    # run msms in with radius 1.0 (if fail, will increase to be 1.1)
    os.system(msmsdir + "msms.x86_64Linux2.2.6.1 -if p_sa.xyzr  -af p_sa.area -probe_radius 1.0 -surface ases > log1.tmp 2>&1")
    os.system(msmsdir + "msms.x86_64Linux2.2.6.1 -if l_sa.xyzr  -af l_sa.area -probe_radius 1.0 -surface ases > log2.tmp 2>&1")
    os.system(msmsdir + "msms.x86_64Linux2.2.6.1 -if pl_sa.xyzr  -af pl_sa.area -probe_radius 1.0 -surface ases > log3.tmp 2>&1")
    if (os.path.isfile('p_sa.area') and os.path.isfile('l_sa.area') and os.path.isfile('pl_sa.area')) == False:
        os.system(msmsdir + "msms.x86_64Linux2.2.6.1 -if p_sa.xyzr  -af p_sa.area -probe_radius 1.1 -surface ases > log1.tmp 2>&1")
        os.system(msmsdir + "msms.x86_64Linux2.2.6.1 -if l_sa.xyzr  -af l_sa.area -probe_radius 1.1 -surface ases > log2.tmp 2>&1")
        os.system(msmsdir + "msms.x86_64Linux2.2.6.1-if pl_sa.xyzr  -af pl_sa.area -probe_radius 1.1 -surface ases > log3.tmp 2>&1")
        print('1.1')


    # read surface area to df2
    df2 = {}
    tmp1 = np.genfromtxt('p_sa.area', skip_header=1)[:,2]
    tmp2 = np.genfromtxt('l_sa.area', skip_header=1)[:,2]
    tmp3 = np.genfromtxt('pl_sa.area', skip_header=1)[:,2]
    df2[2] = tmp2
    df2[3] = tmp3[len(tmp1):]
    df2 = pd.DataFrame(df2)
    df = pd.concat([df1, df2], axis=1)
    df.columns = ['atm','pharma','l','c']
    os.chdir('../')
    os.system('rm -rf tmp_dl')
    os.chdir(olddir)
    return df


def featureSASA(datadir, inprot, inlig, write=False):
    pharmatype = ['P', 'N', 'DA', 'D', 'A', 'AR', 'H', 'PL', 'HA']
    outdict = {i:0 for i in pharmatype}
    df = runMSMS(inprot, inlig, datadir)
    df["d"] = (df["l"] - df["c"]).clip(0,None)
    dfg =  df.groupby("pharma")["d"].sum()
    dfgdict =  dfg.to_dict()
    for i in dfgdict:
        outdict[i] = dfgdict[i]
    sasalist = []
    for i in pharmatype:
        sasalist.append(outdict[i])

    sasalist.append(sum(sasalist))

    if write:
        f = open("sasa.dat", "w")
        f.write(" ".join([str(np.round(i,2)) for i in sasalist]) + "\n")
        f.close()

    return df, sasalist


class sasa_dl:
    def __init__(self, datadir, prot, lig):
        self.datadir = datadir
        self.prot = prot
        self.lig = lig

        self.rawdata, self.sasalist = featureSASA(self.datadir, self.prot, self.lig)

        self.sasaTotal = self.sasalist[-1]
        self.sasaFeatures = self.sasalist[0:-1]

    def info(self):
        featureInfo = ['P', 'N', 'DA', 'D', 'A', 'AR', 'H', 'PL', 'HA']
        return featureInfo

