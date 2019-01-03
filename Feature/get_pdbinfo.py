#! /usr/bin/python
import numpy as np

def isAtom(line):
	if line[0:6] == "ATOM  " or line[0:6] == "HETATM":
		return  True
	else:
		return  False

def isPAtom(line):
    polar_atoms = ["N","O","S"]

    if isAtom(line) and atmn(line).strip()[0] in polar_atoms:
        return True
    else:
        return False

def isIons(line):
    ions = ["MG","MN","CA","ZN","FE","HG","NI","MG","CO","CD","FE2","NA","K","CU","SR","CS","LI"]

    if isAtom(line) and resn(line).strip() in ions:
        return True
    else:
        return False
	
def atmi(line):
    '''
    atom index

    '''

    return line[6:11]

def atmn(line):
    '''
    atom name

    '''

    return line[12:16]

def resn(line):
    '''
    residue name

    '''

    return line[17:20]

def chid(line):
    '''
    chain ID

    '''

    return line[21:22]

def resi(line):
    '''
    residue index

    '''

    return line[22:26]

## combine chid and resi in case resi is larger than 10000

def seqi(line):
    '''
    sequence index

    '''

    return line[21:26]


def coord(line):
    '''
    coordinates

    '''

    crd = [float(line[30 + 8 * i : 38 + 8 * i]) for i in range(3) ]
    return crd

def isHydrogen(line):

    if line[12] == 'H' or line[13] == 'H':
        return 1
    else:
        return 0

def isWater(line):
    if resn(line) == "WAT" or resn(line) == "HOH":
        return True
    else:
        return False



class pdbinfo:

    def __init__(self, name = None, file = None, lines = None):
        self.name = name
        if file != None:
            self.file = file
            self.lines = [line for line in open(self.file)]
        else:
            self.lines = lines

    def getAtoms(self):

        Atoms = [line for line in self.lines if isAtom(line)]

        return Atoms

    def getPolarAtoms(self):
        
        Atoms = [line for line in self.lines if isPAtom(line)]

        return Atoms
    
    def getIons(self):

        Ions = [line for line in self.lines if isIons(line)]

        return Ions

    
    def getProteinWaters(self):

        Waters = [line for line in self.lines if isWater(line)]
        Protein = [ line for line in self.lines if line not in Waters]

        return Protein, Waters

    def getCoords(self):

        Coords = np.array([coord(atom) for atom in self.lines])

        return Coords





    


        
    
