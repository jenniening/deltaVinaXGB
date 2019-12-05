import os
import sys
import get_pdbinfo as pdbinfo
import chimera
from chimera import runCommand as rc
import fileinput
import copy

##### Part1 transfer Hg to Flag or Flag to Hg


def Hg2toFlag(filename):
    """
	change Hg2+ ion atom name to be FLAG
	"""


    f = open("tmp","w")
    for lines in fileinput.input(filename):
        if pdbinfo.isAtom(lines) == 1 and \
            pdbinfo.resn(lines).split()[0] == 'HG':
            lines = lines[0:12] + "FLAG" + lines[16::]
        f.write(lines)
    fileinput.close()
    f.close()
    os.system("mv tmp " + filename)


def FlagtoHg2(filename):
    """
	change FLAG to be HG
	"""
    cmd1 = """sed 's/FLAG/\ HG\ /g' """ + filename + " > tmp"
    cmd2 = "mv tmp " + filename
    os.system(cmd1)
    os.system(cmd2)

##### Part2 use pdb4amber to renumber the residue and atom, and remove the hydrogen

def runPDB4Amber_1(fn):
    pdborg = fn + "_protein.pdb"
    pdbin = fn + "_protein_proc0.pdb"
    os.system("cp " + pdborg + " " + pdbin)
    Hg2toFlag(pdbin)
    pdbout = fn + "_protein_proc1.pdb"
    cmd = "pdb4amber -i " + pdbin + " -o " + pdbout + " --nohyd --dry 2> pdb4amber1.log"
    os.system(cmd)
    FlagtoHg2(pdbin)
    FlagtoHg2(pdbout)

##### Part3 consider the PCA molecule and remove this molecule

def PCA(fn,residue_PCA):
    """
    get the name of files which have PCA

    """
    inputfile =open(fn + "_protein.pdb")
    for line in inputfile:
        if pdbinfo.isAtom(line) == 1 and pdbinfo.resn(line).split()[0] == "PCA":
            residue_PCA.append(fn)

def rmPCA(fn):
    """
    remove PCA from these files

    """
    inputfile = open(fn + "_protein.pdb")
    outputfile = open(fn + "_protein_tmp.pdb","w")
    for line in inputfile:
        if pdbinfo.isAtom(line) == 1 and pdbinfo.resn(line).split()[0] == "PCA":
            continue
        else:
            outputfile.write(line)
    outputfile.close()
    os.system("mv " + fn + "_protein_tmp.pdb " + fn + "_protein.pdb" )
    os.system("rm " + fn + "_tmp.pdb")

##### Part4 keep chains in 12 A of ligand : get fn_protein_proc2.pdb

def getPocketChain(fn):
    """
    get pocket chain

    """
    protpdb = fn + "_protein_proc1.pdb"
    ligpdb = fn + "_ligand.mol2"
    outpdb = fn + "_protein_pocket_12.pdb"
    rc("open " + protpdb)
    rc("open " + ligpdb)
    rc("sel #1 za < 12 & ~#1")
    rc("del ~sel")
    rc("write selected format pdb #0 " + outpdb )
    rc("del #0")

def rmChain(fn):
    """
    remove chains outside 12 A

    """
    pktpdb  = fn + "_protein_pocket_12.pdb"
    protpdb = fn + "_protein_proc1.pdb"
    outpdb  = fn + "_protein_proc2.pdb"
    chidlist = []
    for lines in fileinput.input(pktpdb):
        if pdbinfo.isAtom(lines) == 1:
            if pdbinfo.chid(lines) in chidlist:
                continue
            else:
                chidlist.append(pdbinfo.chid(lines))
    fileinput.close()
    f = open(outpdb,"w")
    for lines in fileinput.input(protpdb):
        if pdbinfo.isAtom(lines) == 1 and pdbinfo.chid(lines) not in chidlist:
            continue
        else:
            f.write(lines)
    fileinput.close()
    f.close()


##### Part5 keep ions in 3 A of protein : get fn_protein_proc3.pdb

def fixFe2(fn):
    """
	Change residue FE2 atom name to be FE2

	"""
    pdbin  = fn + "_protein_proc2.pdb"
    pdbout  = fn + "_protein_proc2_fixFe2.pdb"
    f = open(pdbout, "w")
    for lines in fileinput.input(pdbin):
        if pdbinfo.isAtom(lines) == 1 and pdbinfo.resn(lines) == 'FE2':
            lines = lines[0:12]  +  ' FE2' + lines[16:]
        f.write(lines)
    f.close()
    fileinput.close()
    os.system("mv " + pdbout + " " + pdbin)

def splitProtIons(fn):
    """
	split protein and ions by pdb4amber with option -p
	xxxx_protein_proc2_noion.pdb             <--- protein
	xxxx_protein_proc2_noion_nonprot.pdb     <--- ions

	"""
    pdbin  = fn + "_protein_proc2.pdb"
    Hg2toFlag(pdbin)
    pdbout = fn + "_protein_proc2_noion.pdb"
    cmd = "pdb4amber -i " + pdbin + " -o " + pdbout + " -p 2> pdb4amber21.log"
    os.system(cmd)
    FlagtoHg2(pdbin) # convert Flag back to hg in input pdb

def getCloseIons(fn):
    """
    only keep the close ions 3 angstrom away from protein (Not protein and ligand)
    fn_protein_proc2_ions <--- protein with ions in 3 A

	"""
    inpdb1 =  fn + "_protein_proc2_noion.pdb"
    inpdb2 =  fn + "_protein_proc2_noion_nonprot.pdb"
    FlagtoHg2(inpdb2)
    outpdb  =  fn + "_protein_proc2_ions.pdb"
    if os.stat(inpdb2)[6] > 0:
        rc("open " + inpdb1)
        rc("open " + inpdb2)
        rc("sel #1 & #0 z > 3")
        rc("del sel")
        rc("combine #0,1 modelID 2")
        rc("del #0,1")
        rc("write format pdb #2 " + outpdb )
        rc("del #2")
    else:
        cmd = "cp " + inpdb1 + " " + outpdb
        os.system(cmd)

def runPDB4Amber_2(fn):
    """
    renumber the residue and atom, and remove the hydrogen

	"""
    pdbin  = fn + "_protein_proc2_ions.pdb"
    pdbout = fn + "_protein_proc3.pdb"
    Hg2toFlag(pdbin)
    cmd = "pdb4amber -i " + pdbin + " -o " + pdbout + " --nohyd 2> pdb4amber22.log"
    os.system(cmd)
    FlagtoHg2(pdbin)
    FlagtoHg2(pdbout)
    FlagtoHg2(fn + "_protein_proc3_nonprot.pdb")

def checkIons(fn):
    """
    check the information of ions
    fn_protein_proc3_nonprot <--- ions

	"""
    inpdb = fn + "_protein_proc3_nonprot.pdb"
    ions = []
    if os.stat(inpdb)[6] > 0:
        for lines in fileinput.input(inpdb):
            atm = pdbinfo.atmn(lines)
            res = pdbinfo.resn(lines)
            ions.append(atm + " " + res)
    return ions


##### Part6 fix SE MSE : get fn_protein_proc_se.pdb

def fixMSE(fn):
    """
    fix SE MSE

	"""
    propdb = fn + "_protein_proc3.pdb"
    outpdb = fn + "_protein_proc_se.pdb"
    f = open(outpdb,"w")
    for lines in fileinput.input(propdb):
        if lines[0:4]  == "HETA" and lines[13:20] == "SE  MSE":
            f.write(lines[0:13] + 'S   MSE' + lines[20:76] + ' S\n')
        else:
            f.write(lines)
    f.close()

##### Part7 fix gap: get fn_protein_proc7.pdb

def checkGap(fn):
    """
	Check the gap information from the pdb4amber22.log file

	"""
    gap = []
    for lines in fileinput.input("pdb4amber22.log"):
        if lines[0:3] == "gap":
            line = lines.split()
            #id1 = line[-3].split("_")[1]
            #id2 = line[-1].split("_")[1]
            ### change when using in prepare pair-test ###
            id1 = line[6]
            id2 = line[-1]
            id1 = "%4d" %(int(id1))
            id2 = "%4d" %(int(id2))
            gap.append(id1)
            gap.append(id2)
    gapb = copy.copy(gap)
    gapinch = [] # for gap within one chain
    fileinput.close()
    inpdb = fn + "_protein_proc_se.pdb"
    if len(gap) > 0:
        gapchid = [ i for i in range(len(gap))]
        for lines in fileinput.input(inpdb):
            if pdbinfo.isAtom(lines) == 1 and pdbinfo.resi(lines) in gap:
                tmpid = gap.index(pdbinfo.resi(lines))
                gapchid[tmpid] = pdbinfo.chid(lines)
                gap[tmpid] = 0
        fileinput.close()
        for i in range(len(gapb)/2):
            if gapchid[2*i] == gapchid[2*i+1] and int(gapb[2*i+1]) - int(gapb[2*i]) == 1: # the residues form the gap are close residues in same chain
                gapinch.append(gapb[2*i])
                gapinch.append(gapb[2*i+1])
    return gapinch

def checkCloseGap(gapinch):
    """
    Remove the residues have gaps at two ends ( the residue in two gaps)

    """
    delres = []
    for i in range(len(gapinch)/2-1):
        if gapinch[2*i+1] == gapinch[2*i+2]:
            delres.append(gapinch[2*i+1])
    gapinchdel = copy.copy(gapinch)
    for i in delres:
        gapinchdel.remove(i)
        gapinchdel.remove(i)
    return delres, gapinchdel

def checkGap2(fn):
    gap = []
    inpdb = fn + "_protein_proc7.pdb"
    for lines in fileinput.input("pdb4amber3.log"):
        if lines[0:3] == "gap":
            line = lines.split()
            #id1 = line[-3].split("_")[1]
            #id2 = line[-1].split("_")[1]
            ### change when using in prepare pair-test ###
            id1 = line[6]
            id2 = line[-1]
            id1 = "%4d" %(int(id1))
            id2 = "%4d" %(int(id2) )
            gap.append(id1)
            gap.append(id2)
    gapb = copy.copy(gap)
    gapinch = []
    fileinput.close()
    if len(gap) > 0:
        gapchid = [ i for i in range(len(gap))]
        for lines in fileinput.input(inpdb):
            if pdbinfo.isAtom(lines) == 1 and pdbinfo.resi(lines) in gap:
                tmpid = gap.index(pdbinfo.resi(lines))
                gapchid[tmpid] = pdbinfo.chid(lines)
                gap[tmpid] = 0
        fileinput.close()
        for i in range(len(gapb)/2):
            if gapchid[2*i] == gapchid[2*i+1] and int(gapb[2*i+1]) - int(gapb[2*i]) > 1:
                gapinch.append(gapb[2*i])
                gapinch.append(gapb[2*i+1])
    return gapinch

def addTER(fn, gapinch):
    inpdb = fn + "_protein_proc_se.pdb"
    outpdb = fn + "_protein_proc4.pdb"
    halfgap = [ gapinch[2*i+1] for i in range(len(gapinch)/2) ]
    f = open(outpdb, "w")
    for lines in fileinput.input(inpdb):
        if pdbinfo.isAtom(lines) == 1 and pdbinfo.resi(lines) in halfgap:
            f.write("TER \n")
            halfgap.remove(pdbinfo.resi(lines))
        f.write(lines)
    fileinput.close()
    f.close()

def addGly(fn, gapinch, delres):
    inpdb = fn + "_protein_proc4.pdb"
    outpdb = fn + "_protein_proc5.pdb"
    rc("open " + inpdb)
    gapid = [int(i) for i in gapinch]
    for i in gapid:
        cmd = "addaa gly," + str(i) + "A :" + str(i)
        rc(cmd)
    if len(delres) > 0:
        for i in delres:
            cmd = "del :" + str(int(i))
            rc(cmd)
    rc("write format pdb #0 " + outpdb)
    rc("del #0")


def changeGly(fn, gapinch):
    """
    Change Gly to be ACE or NHE
    """
    inpdb = fn + "_protein_proc5.pdb"
    outpdb = fn + "_protein_proc6.pdb"
    nhe = [gapinch[2*i] + 'A' for i in range(len(gapinch)/2)  ]
    ace = [gapinch[2*i+1] + 'A' for i in range(len(gapinch)/2)  ]
    f = open(outpdb, "w")
    for lines in fileinput.input(inpdb):
        if pdbinfo.isAtom(lines) == 1 and lines[22:27] in nhe:
            if pdbinfo.atmn(lines) == ' N  ':
                lines = lines[:17] +  'NHE' + lines[20:]
            else:
                continue
        elif pdbinfo.isAtom(lines) == 1 and lines[22:27] in ace:
            if pdbinfo.atmn(lines) == ' C  ' or pdbinfo.atmn(lines) == ' O  ':
                lines = lines[:17] +  'ACE' + lines[20:]
            elif pdbinfo.atmn(lines) == ' CA ':
                lines = lines[:17] +  'ACE' + lines[20:]
                lines = lines[:12] +  ' CH3' + lines[16:]
            else:
                continue
        f.write(lines)
    fileinput.close()
    f.close()
    Hg2toFlag(outpdb)
    outpdb2 = fn + "_protein_proc7.pdb"
    cmd = "pdb4amber -i " + outpdb + " -o " + outpdb2 + " --nohyd 2> pdb4amber3.log"
    os.system(cmd)
    FlagtoHg2(outpdb)
    FlagtoHg2(outpdb2)
    gapinch = checkGap2(fn)
    return gapinch

##### Part7 fix missing heavy atom: get fn_protein_proc8.pdb, fn_protein_proc8_water.pdb --> water molecules

def runtleap(inpdb, outpdb1, outpdb2, option):
    f=open("tleap.in","w")
    f.write("""source leaprc.protein.ff14SB
            loadamberparams frcmod.ions1lm_1264_tip3p
            mol = loadpdb %s
            savepdb mol %s
            quit
            """ %(inpdb, outpdb1))
    f.close()
    cmd = "tleap -s -f  tleap.in > tleap.out"
    os.system(cmd)
    Hg2toFlag(outpdb1)
    cmd = "pdb4amber -i " + outpdb1 + " -o " + outpdb2 + " " + option + "2> pdb4amber4.log"
    os.system(cmd)
    FlagtoHg2(outpdb1)
    FlagtoHg2(outpdb2)

##### Part8 Run propka 3.1

def getComplex(fn):
    ## change HIE to HIS after tleap process and remove ACE and NHE
    ## Original change CYX to be CYS, not necessary for propka3.1
    inpdb = fn + "_protein_proc8.pdb"
    outpdb1 = fn + "_protein_proc8_noCYX.pdb"
    cmd1 = """sed 's/HIE/HIS/g;s/HID/HIS/g;s/HIP/HIS/g;/ACE/d;/NHE/d' """ + inpdb + " > " + outpdb1
    os.system(cmd1)
    ## read in ligand mol2, fix the atom naming by convert to ac
    ## then convert back to mol2
    ## atom type should be sybyl since Chimera cannot recognize gaff
    ## (Dont use PDB, it will make mistake of element)
    inmol2 = fn + "_ligand.mol2"
    outmol2 = fn + "_ligand_fixed.mol2"
    outac = fn + "_ligand_fixed.ac"
    cmd =  "antechamber -fi mol2 -fo ac -i " + inmol2 + " -o " + outac + " -at sybyl -dr no"
    os.system(cmd)
    cmd =  "antechamber -fi ac -fo mol2 -i " + outac + " -o " + outmol2 + " -at sybyl -dr no"
    os.system(cmd)
    ## get complex of protein and ligand##
    outpdb2 = fn + "_complex_proc9.pdb"
    rc("open " + outpdb1)
    rc("open " + outmol2)
    rc("combine #0,1 modelID #2")
    rc("write format pdb #2 " + outpdb2)
    rc("close all")

def runPropka(fn):
    inpdb = fn + "_complex_proc9.pdb"
    cmd = propka31 + " -q " + inpdb + " > propka31.log"
    os.system(cmd)


#####  Part 9 Correct protein protonation state 
##### Read the propka31 information from .pka file, decide the protonation state of HIS/HIP, ASP/ASH, GLU/GLH, LYS/LYH, CYS/CYM
##### Run pdb2pqr to get protonation state of HIS, ASP, GLU, LYS
##### Input: xxxx_protein_proc8.pdb
##### Output: xxxx_protein_proc8_noCYX.pdb, xxxx_protein_proc8_noCYX.pqr, xxxx_protein_prep.pdb

def _checkProtpka(modelpka, calcpka):
    """
    status = flag1
    if modelpka and calcpka on same side of 7: 0
    if modelpka and calcpka on two side of 7: 1

    """
    modelpka = float(modelpka)
    calcpka = float(calcpka)
    if modelpka >= 7 and calcpka >= 7:
        status = 0
    elif modelpka >= 7 and calcpka < 7:
        status = 1
    elif modelpka < 7 and calcpka >= 7:
        status = 1
    elif modelpka < 7 and calcpka < 7:
        status = 0
    return status

def readpKa(fn):
    infn = fn + "_complex_proc9.pka"
    res = {}
    resdict = {"   ASP":"ASH", "   GLU":"GLH", "   HIS":"HIP", "   LYS":"LYN", "   CYS":"CYM"}
    for lines in fileinput.input(infn):
        if lines[0:6] in resdict.keys():
            line = lines[6::].split()
            status = _checkProtpka(line[3], line[2])
            if status == 1:
                res[line[0]]= resdict[lines[0:6]]
    fileinput.close()
    return res

def assignPDB1(fn, res):

    chid = res.keys()
    inpdb = fn + "_protein_proc8_noCYX.pdb"
    outpdb = fn + "_protein_proc8_assign1.pdb"
    f = open(outpdb, "w")
    for lines in fileinput.input(inpdb):
        if pdbinfo.isAtom(lines) == 1:
            resid = pdbinfo.resi(lines).split()[0]
            if resid in chid:
                lines = lines[:17] +  res[resid] + lines[20:]
            f.write(lines)
    fileinput.close()
    f.close()

def runPDB2PQR(fn):
    inpdb = fn + "_protein_proc8_assign1.pdb"
    outpdb = fn + "_protein_proc8_assign1.pqr"
    cmd = pdb2pqr + " --ff=amber  -v " + inpdb + " " + outpdb + " > pdb2pqr.log "
    os.system(cmd)

def checkProtonState(proton):
    res = {}
    for i in proton.keys():
        if proton[i][0] == 'HIS':
            if len(proton[i][1]) == 2:
                res[i] = 'HIP'
            elif proton[i][1][0] == ' HD1':
                res[i] = 'HID'
            elif proton[i][1][0] == ' HE2':
                res[i] = 'HIE'
    return res

def getProton(fn):
    inpdb = fn + "_protein_proc8_assign1.pqr"
    proton = {}
    for lines in fileinput.input(inpdb):
        if pdbinfo.isAtom(lines) == 1:
            resn = pdbinfo.resn(lines)
            if resn == 'HIS':
                resi = pdbinfo.resi(lines).split()[0]
                atmn = pdbinfo.atmn(lines)
                if resi in proton.keys() and atmn in [' HD1', ' HE2']:
                    proton[resi][1].append(atmn)
                elif resi not in proton.keys():
                    proton[resi] = [resn,[]]
                    if atmn in [' HD1', 'HE2']:
                        proton[resi][1].append(atmn)
    res = checkProtonState(proton)
    return res

def changeRes(fn, res1, res2):
    inpdb = fn + "_protein_proc8.pdb"
    outpdb = fn + "_protein_prep.pdb"
    chid1 = res1.keys()
    chid2 = res2.keys()
    f = open(outpdb, "w")
    for lines in fileinput.input(inpdb):
        if pdbinfo.isAtom(lines) == 1:
            resid = pdbinfo.resi(lines).split()[0]
            if resid in chid1:
                lines = lines[:17] +  res1[resid] + lines[20:]
            elif resid in chid2:
                lines = lines[:17] +  res2[resid] + lines[20:]
            f.write(lines)
    fileinput.close()
    f.close()




class general_clean_pro:
    def __init__(self, olddir, datadir, pllist, residue_PCA=None):
        """
        clean protein
        
        :param olddir: current directory
        :param datadir: directory for structure file
        :param pllist: the list of all pdb index 
        :residue_PCA: the list of pdb indexes that have PCA, defaults to None

        """

        self._olddir = olddir
        self._datadir = datadir
        self._pllist = pllist
        self._residue_PCA = residue_PCA
        self._out = os.path.join(datadir, "gap_notfix_pep.txt")
        self._out2 = os.path.join(datadir, "gap_check_pep.txt")
        self._out_proton = os.path.join(datadir, "protonation_add.txt")
        self._out_problem = os.path.join(datadir, "prep_pro_problem.txt")
        os.system("""echo 'PDBID # gaps' > """ + self._out)
        os.system("""echo 'PDBID # gaps' > """ + self._out2)

    @property
    def datadir(self):
        return self._datadir

    @property
    def olddir(self):
        return self._olddir

    @property
    def pllist(self):
        return self._pllist

    @property
    def out(self):
        return self._out

    @property
    def out2(self):
        return self._out2
    
    @property
    def out_proton(self):
        return self._out_proton
    
    @property
    def residue_PCA(self):
        if self._residue_PCA == None:
            residue_PCA = []
            for fn in self.pllist:
                os.chdir(os.path.join(self.datadir,fn))
                PCA(fn, residue_PCA)
            self._residue_PCA = set(residue_PCA)
        return self._residue_PCA
    
    def removePCA(self):
        for fn in self.residue_PCA:
            os.chdir(os.path.join(self.datadir,fn))
            rmPCA(fn)
        
    def cleanProp_1(self):
        for fn in self.pllist:
            os.chdir(os.path.join(self.datadir, fn))
            os.system("rm pdb4amber*")
            runPDB4Amber_1(fn)
            getPocketChain(fn)
            rmChain(fn)
            fixFe2(fn)
            splitProtIons(fn)
            getCloseIons(fn)
            runPDB4Amber_2(fn)
            fixMSE(fn)
            ions = checkIons(fn)
            numions = len(ions)
            out = fn +'_protein_ions_in3_pep.pdb'
            if out in os.listdir(os.path.join(self.datadir, fn)):
                os.system("rm " + out)
            if  numions > 0:
                cmd = """echo '""" + fn + ' ' + str(numions) + ' ' + ''.join(ions) +  """'  >> """ + out
                os.system(cmd)

    def cleanProp_2(self):
        for fn in self.pllist:
            os.chdir(os.path.join(self.datadir,fn))
            os.system("echo " + fn)
            gapinch = checkGap(fn)
            numgap = len(gapinch)/2
            if numgap > 0:
                cmd = """echo '""" + fn + ' ' + str(numgap) + ' ' + ''.join(gapinch) +  """'  >> """ + self.out
                os.system(cmd)
                if numgap > 1:
                    delres, gapinchdel = checkCloseGap(gapinch)
                else:
                    delres = []
                    gapinchdel = copy.copy(gapinch)
                addTER(fn, gapinch)
                addGly(fn, gapinchdel, delres)
                gapinch1 = changeGly(fn, gapinchdel)
                if len(gapinch1) > 0:
                    cmd = """echo '""" + fn + ' ' + ''.join(gapinch1) +  """'  >> """ + self.out2
                    os.system(cmd)
            else:
                cmd = "cp " +  fn + "_protein_proc3.pdb  " + fn + "_protein_proc7.pdb"
                os.system(cmd)
            inpdb = fn + "_protein_proc7.pdb"
            outpdb1 = fn + "_protein_proc7_addh.pdb"
            outpdb2 = fn + "_protein_proc8.pdb"
            option = "--nohyd --dry "
            runtleap(inpdb, outpdb1, outpdb2, option)
            os.chdir(self.olddir)
    
    def run_proka(self):
        for fn in self.pllist:
            os.chdir(os.path.join(self.datadir,fn))
            os.system("echo " + fn)
            getComplex(fn)
            runPropka(fn)

    def correct_proton(self):
        for fn in self.pllist:
            os.chdir(self.datadir + fn)
            if fn + "_complex_proc9.pka" in os.listdir(os.path.join(self.datadir, fn)):
                res1 = readpKa(fn)
                assignPDB1(fn, res1)
                runPDB2PQR(fn)
                res2 = getProton(fn)
                tmp = ''
                for i in [res1, res2]:
                    if len(i) > 0:
                        for j in i.keys():
                            tmp = tmp  + " " + j + " " + i[j]
                cmd = """echo '""" + fn + ' ' + tmp +  """'  >> """ + self.out_proton
                os.system(cmd)
                changeRes(fn, res1, res2)
                inpdb = fn + "_protein_prep.pdb"
                outpdb1 = fn + "_protein_prep_proton.pdb"
                outpdb2 = fn + "_protein_prep_final.pdb"
                option = ""
                runtleap(inpdb, outpdb1, outpdb2, option)
            else:
                print(fn, "Wrong propka")


if __name__ == "__main__":
    ### directory for propka31 and pdb2pqr ###
    propka31 = "/Users/jianinglu1/.local/bin/propka31"
    pdb2pqr = "/Users/jianinglu1/Documents/Tools/pdb2pqr-osx-bin64-2.1.0/pdb2pqr"

    datadir = ""
    olddir = os.getcwd()
    pllist = ["2al5"]
    propclean = general_clean_pro(olddir, datadir, pllist)
    propclean.cleanProp_1()
    propclean.cleanProp_2()
    propclean.run_proka()
    propclean.correct_proton()




