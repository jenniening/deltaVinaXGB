"""
Pharmacophore type assignment for protein and ligand.
"""

__author__ = "Jianing Lu"
__copyright__ = "Copyright 2017, NYU"
__license__ = ""

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import os, sys
import pybel
import openbabel as ob

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------

class pharma:
    """Pharmacophore Type Assignment for Protein and Ligand
    
    This class is to assign the pharmacophore type 
    
    Parameters
    ----------
    fn : str
        File name of protein or ligand
    
    References
    ----------
    .. [1] Our paperReference of Dock
    .. [2] Reference of Dock

    """
    def __init__(self, fn):
        self.fn = fn
        
        # Empty list for AtomIdx and Empty dict for pharma
        self.AtomIdx = []
        self.AtomPharma = {}


    def assign(self, write = False, outfn = 'tmp.pdb'):
        """Assign pharmacophore type.
        
        Details description of pharmacophore assign are commented in code.
        
        Parameters
        ----------
        write : logic 
            Control to write the filtered PDB after pharmacophore assignment
            Remove atoms not in element list
        outfn : str
            Output PDB file name with default 'tmp.pdb'
            
        Returns
        ---------
        AtomIdx : list
            Atom index in the structure
        AtomPharma : dict
            Dict with atom index as key and pharmacophore type as value
        
        """
        # Nine element will be used in the study
        elementint = [6, 7, 8, 9, 15, 16, 17, 35, 53]

        # supress the logging information 
        pybel.ob.obErrorLog.StopLogging()
        
        # table of convert OB internal atom type to Sybyl
        ttab = pybel.ob.OBTypeTable()
        ttab.SetFromType("INT")
        ttab.SetToType("SYB")
        
        # read in molecule
        __, ft = os.path.splitext(self.fn)
        mol = pybel.readfile(ft[1:], self.fn).__next__()

        # convert the atom type from internal to sybyl
        for atom in mol.atoms:
            # AtomIdx.append(atom.idx)
            # convert the atom type and make it upper case
            at = ttab.Translate(atom.OBAtom.GetType())
            at = at.upper()
            atom.OBAtom.SetType(at)
        
        # assign pharmacophore type
        for atom in mol.atoms:
            # append atom idx to the AtomIdx
            self.AtomIdx.append(atom.idx)
            at = atom.type
            #print at

            # pharma type for element not in C,N,O,P,S,F,Cl,Br,I
            if atom.atomicnum not in elementint:
                p = 'NU'
            
            # pharma type for oxygen
            elif at in ['O.3', 'O.2', 'O.CO2']:
                p = 'A'
                # nbr of oxygen to be one atom or error
                nbrs = [i for i in ob.OBAtomAtomIter(atom.OBAtom)]
                if len(nbrs) == 2:
                    for nbr in ob.OBAtomAtomIter(atom.OBAtom):
                        if nbr.GetAtomicNum() == 1:
                            p = 'DA'
                            
                elif len(nbrs) == 1:
                    nbr = list(nbrs)[0]
                    
                    # nbr is carbon check if it is coo-
                    if nbr.GetAtomicNum() in [6, 15]:
                        c = 0
                        for nbr2 in ob.OBAtomAtomIter(nbr):
                            if nbr2.GetAtomicNum() in [8,16]:
                                if len(list(ob.OBAtomAtomIter(nbr2))) == 1:
                                    c += 1
                        if c >= 2:
                            p = 'N'
                            
                    elif nbr.GetAtomicNum() == 16:
                        c = 0
                        for nbr2 in ob.OBAtomAtomIter(nbr):
                            if nbr2.GetAtomicNum() == 8:
                                if len(list(ob.OBAtomAtomIter(nbr2))) == 1:
                                    c += 1
                        if c >= 3:
                            p = 'N'             
            
            # pharma type for nitrogen              
            elif at == 'N.4':
                p = 'P'
            elif at == 'N.3':
                p = 'A'
                for nbr in ob.OBAtomAtomIter(atom.OBAtom):
                    if nbr.GetAtomicNum() == 1:
                        p = 'DA'
                        break
            elif at == 'N.2':
                p = 'A'
                nbrs = [i for i in ob.OBAtomAtomIter(atom.OBAtom)]
                if len(nbrs) == 3:
                    p = 'P'
                else:
                    for nbr in ob.OBAtomAtomIter(atom.OBAtom):
                        if nbr.GetAtomicNum() == 1:
                            p = 'DA'
            elif at == 'N.1':
                p = 'A'
            elif at == 'N.AR':
                p = 'AR'
                nbrs = [i for i in ob.OBAtomAtomIter(atom.OBAtom)]
                if len(nbrs) == 3:
                    for nbr in ob.OBAtomAtomIter(atom.OBAtom):
                        if nbr.GetAtomicNum() == 1:
                            p = 'D'
                elif len(nbrs) == 2:
                    p = 'A'
            elif at == 'N.AM':
                p = 'PL'
                for nbr in ob.OBAtomAtomIter(atom.OBAtom):
                    if nbr.GetAtomicNum() == 1:
                        p = 'D'           
            elif at == 'N.PL3':
                p = 'A'
                for nbr in ob.OBAtomAtomIter(atom.OBAtom):
                    #print atom.idx, at, nbr.GetAtomicNum() ,atom.OBAtom.GetBond(nbr).GetBondOrder()
                    if nbr.GetType() == 'C.CAT':
                        p = 'P'
                        break  
                    elif nbr.GetAtomicNum() == 1:
                        p = 'DA'      
                        
            # pharma type for sulfur    
            elif at in ['S.3', 'S.2', 'S.O', 'S.O2']:
                p = 'PL'
                nbrs = [i for i in ob.OBAtomAtomIter(atom.OBAtom)]
                if len(nbrs) == 1:
                    p = 'A'
                    nbr = nbrs[0]           
                    if nbr.GetAtomicNum() == 6:
                        nbrs2 = [i for i in ob.OBAtomAtomIter(nbr)]
                        if len(nbrs2) == 4:
                            p = 'N'
                        elif len(nbrs2) == 3:
                            c = 0
                            for nbr3 in ob.OBAtomAtomIter(nbr):
                                if nbr3.GetAtomicNum() in [8,16]:
                                    if len(list(ob.OBAtomAtomIter(nbr3))) == 1:
                                        c += 1
                            if c >= 2:
                                p = 'N'        
                elif len(nbrs) == 2:
                    p = 'A'
                    for nbr in ob.OBAtomAtomIter(atom.OBAtom):
                        if nbr.GetAtomicNum() == 1:
                            p = 'DA'    
                    
            # pharma type for carbon              
            elif at == 'C.AR':
                p = 'AR'
            elif at in ['C.1','C.2', 'C.3', 'C.CAT']:
                p = 'H'
                for nbr in ob.OBAtomAtomIter(atom.OBAtom):
                    if nbr.GetAtomicNum() in [7, 8, 9, 15, 16]:
                        p = 'PL'
                        break          
            
            # pharma type for P and Halogen
            elif at == 'P.3':
                p = 'PL'
            elif at in ['F', 'CL', 'BR', 'I']:
                p = 'HA'
            
            # pharma type for general carbon not be assigned
            elif atom.atomicnum == 6:
                p = 'H'
            # pharma type for N,O,F,S,P,Cl,Br,I not be assigned
            elif atom.atomicnum in elementint:
                p = 'PL'
            
            # AtomPharma dict with atomicnum, pharma tyep, and coords
            # the coords is for SADE only
            self.AtomPharma[atom.idx] = [atom.atomicnum, p, atom.coords]
            #print atom.idx, AtomPharma[atom.idx]
            #print atom.idx, atom.type,  p, atom.OBAtom.GetResidue().GetName()

        if write:
            for idx in self.AtomIdx[::-1]:
                if self.AtomPharma[idx][0] not in elementint:
                    mol.OBMol.DeleteAtom(mol.OBMol.GetAtom(idx))
            
            output = pybel.Outputfile("pdb", outfn, overwrite=True)
            output.write(mol)
            output.close()
            
        return self.AtomIdx, self.AtomPharma


    def writePharma(self):
        """Write PDB with pharmacophore type with extension .pharma.pdb
        
        TODO
        ----------
        * Support write .pharm.pdb with mol2 input
        
        """
        AtomIdx, AtomPharma = pharmaAssign(fn)
        
        # new pdb file
        ftmp, fext = os.path.splitext(self.fn)
        
        if fext[1:].lower() != 'pdb':
            print("Sorry, currently only support input PDB to write.")
            sys.exit()
        
        outfn = ftmp + ".pharma.pdb"

        f = open(outfn, "w")
        with open(fn) as g:
            for lines in g:
                if lines[0:6] in ['ATOM  ', 'HETATM']:
                    atmi = int(lines[6:11])
                    #reschain = lines[17:26]
                    lines = lines[0:76] + "%2s" + lines[78:]
                    f.write(lines %(self.AtomPharma[atmi]))
                else:
                    f.write(lines)
        f.close()



def test():

    print("""Test: protein pocket pharmacophore for write implementation """)
    
    datadir = "/Users/chengwang/Dropbox/ws/scoref/v13/test/"
    pdblist = ['1a8i']
    for pdb in pdblist:
        fn = datadir + pdb + '/1a8i_protein.pdb'
        AtomIdx, AtomPharma = pharmaAssign(fn, write=True)
        fn = datadir + pdb + '/1a8i_decoys_native.mol2'
        AtomIdx, AtomPharma = pharmaAssign(fn, write=True)     

    print("""End of test4""")



if __name__ == '__main__':
    
    test()
