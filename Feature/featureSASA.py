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
    """Assign pharmaphore type to each atom and calculate SASA by MSMS

    Details can be found in comments

    Parameters
    ----------
    inprot : str
        Protein structure input
    inlig : str
        Ligand structure input

    Returns
    ----------
    df

    """

    # create tmp folder for all intermediate files
    olddir = os.getcwd()
    os.chdir(MSMSDIR)
    os.system('mkdir tmp')

    # Process Input file to be p.pdb and l.pdb
    # convert protein file to PDB if not and remove hetatm card
    ppdb = 'tmp/p.pdb'
    __, intype = os.path.splitext(inprot)
    if intype[1:].lower() != 'pdb':
        prot = pybel.readfile(intype[1:], inprot).__next__()
        output = pybel.Outputfile("pdb", ppdb, overwrite=True)
        output.write(prot)
        output.close()
    else:
        # change possible HETATM to ATOM in pdb
        os.system("""sed 's/HETATM/ATOM\ \ /g' """ + inprot + " > " + ppdb)

    # convert ligand file to be PDB by openbabel
    lpdb = 'tmp/l.pdb'
    __, intype = os.path.splitext(inlig)
    if intype[1:].lower() != 'pdb':
        lig = pybel.readfile(intype[1:], inlig).__next__()
        output = pybel.Outputfile("pdb", lpdb, overwrite=True)
        output.write(lig)
        output.close()
    else:
        # change possible HETATM to ATOM in pdb
        os.system("""sed 's/HETATM/ATOM\ \ /g' """ + inlig + " > " + lpdb)

    os.chdir('tmp')
    #copy atom typefiel into directory
    os.system("cp " + msmsdir + "atmtypenumbers .")
    # Process p.pdb/l.pdb to be p_sa.pdb/l_sa.pdb after pharma assignment
    # get full atom idx list and pharma
    ppdb2 = 'p_sa.pdb'
    lpdb2 = 'l_sa.pdb'

    pidx, ppharm = pharma('p.pdb').assign(write=True, outfn = ppdb2)
    lidx, lpharm = pharma('l.pdb').assign(write=True, outfn = lpdb2)

    # get subset atom idx which is nine element type
    # This have been done in pharma but still do it again
    elementint = [6, 7, 8, 9, 15, 16, 17, 35, 53]
    psub = [idx for idx in pidx if ppharm[idx][0] in elementint]
    lsub = [idx for idx in lidx if lpharm[idx][0] in elementint]

    # get element number and pharma type and assign to df1
    comp = []
    for idx in psub:
        comp.append(ppharm[idx][0:2])
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
    df2[2] = np.append(tmp1, tmp2)
    df2[3] = tmp3
    df2 = pd.DataFrame(df2)

    df = pd.concat([df1, df2], axis=1)
    df.columns = ['atm','pharma','pl','c']
    os.chdir('../')
    os.system('rm -rf tmp')
    os.chdir(olddir)
    return df


def featureSASA(datadir, inprot, inlig, write=False):
    """Group the SASA by pharmacophore type

    Details can be found in comments

    Parameters
    ----------
    inprot : str
        Protein structure input
    inlig : str
        Ligand structure input

    Returns
    ----------
    sasalist : list [float]

    """

    # nine elements and nine pharma types
    #elemint = [6, 7, 8, 9, 15, 16, 17, 35, 53]
    #elemstr = [str(i) for i in elemint]
    pharmatype = ['P', 'N', 'DA', 'D', 'A', 'AR', 'H', 'PL', 'HA']
    outdict = {i:0 for i in pharmatype}

    # run MSMS
    df = runMSMS(inprot, inlig, datadir)

    ## delta SASA with clip 0 (if value less 0, cut to 0)
    df["d"] = (df["pl"] - df["c"]).clip(0,None)

    # group delta sasa by element and pharma type
    dfg =  df.groupby("pharma")["d"].sum()
    dfgdict =  dfg.to_dict()

    # assign grouped dict to outdict
    for i in dfgdict:
        outdict[i] = dfgdict[i]

    # output list
    sasalist = []
    for i in pharmatype:
        sasalist.append(outdict[i])

    sasalist.append(sum(sasalist))

    if write:
        f = open("sasa.dat", "w")
        f.write(" ".join([str(np.round(i,2)) for i in sasalist]) + "\n")
        f.close()

    return df, sasalist


class sasa:
    """Buried SASA features

    """

    def __init__(self, datadir,prot, lig):
        """Pharmacophore based buried SASA Features

        Parameters
        ----------
        prot : str
            protein structure
        lig : str
            ligand structure

        """
        self.prot = prot
        self.lig = lig
        self.datadir = datadir

        self.rawdata, self.sasalist = featureSASA(self.datadir, self.prot, self.lig)

        self.sasaTotal = self.sasalist[-1]
        self.sasaFeatures = self.sasalist[0:-1]

    def info(self):
        """Feature list"""
        featureInfo = ['P', 'N', 'DA', 'D', 'A', 'AR', 'H', 'PL', 'HA']
        return featureInfo
