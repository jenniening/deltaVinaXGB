##################################################
##### BACE1 cyclic ligand pose fragmentation
##### Author: cyang; 
##### Time: 11/16/2018
##### Requirement: numpy, rdkit
##################################################

import os, sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

def get_spec_N(mol,AtomIdxs):
    ''' 
    Get the N index for N with no less than 3 bonds

    '''

    spec_N_type1 = []
    spec_N_type2 = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetIdx() in AtomIdxs and atom.IsInRing() == True:
            if len(atom.GetBonds()) == 3 and atom.GetFormalCharge() == 0:
                spec_N_type1.append(atom.GetIdx())
                print("There is a spec N type1:" + str(atom.GetIdx()))
            elif len(atom.GetBonds()) == 3 and atom.GetFormalCharge() != 0:
                spec_N_type2.append(atom.GetIdx())
                print("There is a spec N type2:" + str(atom.GetIdx()))
            elif len(atom.GetBonds()) == 4:
                spec_N_type1.append(atom.GetIdx())
                print("There is a spec N type1:" + str(atom.GetIdx()))

            

    return spec_N_type1, spec_N_type2


def check_cutting(rings, max_size_idx, mol):
    '''
    Check the number of sidechains after cutting, more is better

    '''

    max_fgs = 0
    for idx in max_size_idx:
        _,sidechain_connection,_ = GetMacroCyclicFragmentIdx(rings, idx, mol=mol, confId=-1)
        if len(sidechain_connection) > max_fgs:
            max_fgs = len(sidechain_connection)
            max_idx = idx
    return max_idx


def max_ring(mol, rings, cut_type = "only1"):
    '''
    Get the large ring index based on size and cutting 
    
    :param cut_type: how many core used, defaults to "only1"

    '''

    size_infor = [len(i) for i in rings]
    max_size = len(rings[0])
    max_size_idx = [idx for idx, i in enumerate(size_infor) if i == max_size]

    if len(max_size_idx) == 1:
        core_idx = [0]

    elif len(max_size_idx) != 1:
        if cut_type == "only1":
            # use the number of generated sidechains to determin, more is better
            core_idx = [check_cutting(rings,max_size_idx,mol)]

        elif cut_type == "all":
            core_idx = max_size_idx

    return core_idx
    

def determine_intersect(ring_sets):
    Flag = False
    for i in ring_sets:
        for j in ring_sets:
            if i != j:
                if i.intersection(j) or i.issubset(j):
                    Flag = True
    return Flag
                

def get_large_ring(cut_type, mol=None, confId=-1):
    '''
    Get ring element infor and the lager ring index
    
    '''
    rings = sorted([set(i) for i in Chem.GetSymmSSSR(mol)], key=lambda x: len(x), reverse=True)

    newrings = []

    for i in range(0,len(rings)):
        newlist = [k for k in rings[i]]
        for j in range(0, len(rings)):
            if set(newlist).intersection(rings[j]):
                newlist.extend(list(rings[j]))    
        if set(newlist) not in newrings:
            newrings.append(set(newlist))
    newrings = sorted(newrings,key=lambda x: len(x), reverse=True)

    n = 0 
    while determine_intersect(newrings):
        latestrings = []
        for i in newrings:
            newlist = [n for n in i]
            for j in newrings:
                if i == j: 
                    continue
                elif j.issubset(i):
                    continue
                elif i.intersection(j):
                    newlist.extend([m for m in j])
            if set(newlist) not in latestrings:
                latestrings.append(set(newlist))
        newrings = latestrings
 
    
    newrings = sorted(newrings,key=lambda x: len(x), reverse=True)
    
    if len(newrings[0]) <= 10:
        # combine all rings connected by one single bond
        print("The largest ring size: " + str(len(newrings[0])))
    
        latestrings = []
        for i in range(0,len(newrings)):
            newlist = [k for k in newrings[i]]
            for k in newrings[i]:
                atom = mol.GetAtomWithIdx(k)
                for n in atom.GetNeighbors():
                    idx = n.GetIdx()
                    for j in range(i+1, len(newrings)):
                        if idx in newrings[j]:
                            newlist.extend(list(newrings[j]))
            latestrings.append(set(newlist))
        newrings = latestrings
        newrings = sorted(newrings,key=lambda x: len(x), reverse=True)
            
        while determine_intersect(newrings):
            latestrings = []
            for i in newrings:
                newlist = [n for n in i]
                for j in newrings:
                    if i == j: 
                        continue
                    elif j.issubset(i):
                        continue
                    elif i.intersection(j):
                        newlist.extend([m for m in j])
                if set(newlist) not in latestrings:
                    latestrings.append(set(newlist))
            newrings = latestrings


    # get large ring index
    biggest_index_list = max_ring(mol, newrings, cut_type = cut_type )
    #biggest = newrings[biggest_index]

    return newrings, biggest_index_list


def GetMacroCyclicFragmentIdx(newrings, biggest_index, mol=None, confId=-1):
    '''
    Get the AtomIdx of core fragment and the connection information with side chains

    '''
    print("The core size after connection: "  + str(len(newrings[biggest_index])))
    biggest = newrings[biggest_index]
    
    biggest_and_neighbor = biggest.copy()
    sidechain_connection = [] ### belong to the side_chain connection 
    core_connection = [] ### belong to the MacroCyclicFragment connection
    for i in biggest:
        atom = mol.GetAtomWithIdx(i)
        for j in atom.GetNeighbors():
            k = j.GetIdx()
            if k not in biggest_and_neighbor:
                bond = [ bond for bond in atom.GetBonds() if (bond.GetBeginAtomIdx() == k) or (bond.GetEndAtomIdx() == k)][0] 
                if bond.GetBondType() == Chem.BondType.SINGLE and j.GetAtomicNum() not in [1] :
                    sidechain_connection.append(k)
                    core_connection.append(i)
                else:
                    biggest_and_neighbor.add(k)
                
    oldAtomIdxs = biggest_and_neighbor
    return biggest_and_neighbor, sidechain_connection, core_connection



def GetFragmentConf(AtomIdxs, connectIdxs, mol=None, confId=-1):
    '''
    Get the 3D conformation of newly edited fragment
    
    '''

    oldAtomIdxs = AtomIdxs
    #print(AtomIdxs)
    #print(connectIdxs)

    ed_mol   = Chem.rdchem.EditableMol(Chem.rdchem.Mol())
    conf_old = mol.GetConformer(confId)
    conf_new = Chem.rdchem.Conformer()
    mapOldIdxToNewIdx = {}
    cutoff_index = []
    n = 0 
    for oldAtom in mol.GetAtoms():
        oldAtomIdx = oldAtom.GetIdx()
        if oldAtomIdx in oldAtomIdxs:
            newAtom = Chem.rdchem.Atom(oldAtom.GetAtomicNum())
            # correct charge
            newAtom.SetFormalCharge(oldAtom.GetFormalCharge())
            # set hybradization for adding H
            newAtom.SetHybridization(oldAtom.GetHybridization())

            newAtomIdx = ed_mol.AddAtom(newAtom)
            mapOldIdxToNewIdx.update({oldAtomIdx: newAtomIdx})
            point3D = conf_old.GetAtomPosition(oldAtomIdx)
            conf_new.SetAtomPosition(newAtomIdx, point3D)
            if oldAtomIdx in connectIdxs:
                cutoff_index.append(n)
            n += 1
    spec_N_type1, spec_N_type2 = get_spec_N(mol,connectIdxs)
    
    for bond in mol.GetBonds():
        bIdx = bond.GetBeginAtomIdx()
        eIdx = bond.GetEndAtomIdx()
        if bIdx in oldAtomIdxs and eIdx in oldAtomIdxs:
            if bIdx in spec_N_type1 or eIdx in spec_N_type1:
                ed_mol.AddBond(mapOldIdxToNewIdx[bIdx], mapOldIdxToNewIdx[eIdx],order = Chem.BondType.SINGLE)
            else:
                ed_mol.AddBond(mapOldIdxToNewIdx[bIdx], mapOldIdxToNewIdx[eIdx],order = bond.GetBondType())
    
    mol = ed_mol.GetMol()
    #Chem.SanitizeMol()
    mol.AddConformer(conf_new)
    # add H for connection atoms
    mol.UpdatePropertyCache()
    #mol_H = Chem.AddHs(mol,explicitOnly=False,addCoords=True,onlyOnAtoms=cutoff_index)

    return mol

def get_atoms_to_visit(atom, seen_ids):
    '''
    search for all the available neightbors in the side chain
    '''
    neighbor_ids = []
    neighbor_atoms = []
    for bond in atom.GetBonds():
        neighbor_atom = bond.GetOtherAtom(atom)
        neighbor_id = neighbor_atom.GetIdx()
        if neighbor_id not in seen_ids:
            neighbor_ids.append(neighbor_id) 
            neighbor_atoms.append(neighbor_atom)
    
    return neighbor_ids, neighbor_atoms

def find_atom_ids_to_delete(mol, start_atom_id, ignore_atom_id):
    '''
    start_atom_id: connecting atom in the side_chain
    ignore_atom_id: connecting atom in the core
    '''
    seen_ids = {start_atom_id, ignore_atom_id}
    atom_ids = [start_atom_id] # starting with first side-chain atom
    stack = [mol.GetAtomWithIdx(start_atom_id)]
    
    # Use a depth-first search to find the connected atoms
    while stack:
        atom = stack.pop()
        atom_ids_to_visit, atoms_to_visit = get_atoms_to_visit(atom, seen_ids)
        atom_ids.extend(atom_ids_to_visit)
        stack.extend(atoms_to_visit)
        seen_ids.update(atom_ids_to_visit)
    atom_ids = sorted(atom_ids)
    return atom_ids

def WriteFragmentSDF(mol, mol_name):
    writer = Chem.SDWriter(mol_name)
    writer.SetKekulize(False)
    writer.write(mol)
    writer.close()


def WriteFragmentPDB(mol, mol_name):
    Chem.MolToPDBFile(mol,mol_name)



def runMethod(mol, basename,cut_type):
    '''
    Run method 

    '''
    rings = sorted([set(i) for i in Chem.GetSymmSSSR(mol)], key=lambda x: len(x), reverse=True)
    if len(rings) == 0:
        print("No Rings")
        WriteFragmentPDB(mol, '%s_frag%02d.sdf'%(basename.split(".")[0], 1))
    elif cut_type == "only1":
        newrings, biggest_index_list = get_large_ring(cut_type, mol=mol)
        coreIdxs, sidechain_connectIdxs, core_connectIdxs = GetMacroCyclicFragmentIdx(newrings, biggest_index_list[0], mol=mol)
        core = GetFragmentConf(coreIdxs, core_connectIdxs, mol, confId=-1)
        count = 1
        #WriteFragmentSDF(core, '%s_frag%02d.sdf'%(basename.split(".")[0], count)) ## write down core.sdf
        WriteFragmentPDB(core, '%s_frag%02d.pdb'%(basename.split(".")[0], count))
        for m, n in zip(sidechain_connectIdxs, core_connectIdxs):
            sidechainIdx = find_atom_ids_to_delete(mol, m, n)
            sidechain = GetFragmentConf(sidechainIdx, [m], mol)
            count = count + 1
            #WriteFragmentSDF(sidechain, '%s_frag%02d.sdf'%(basename.split(".")[0],count)) ## write down sidechain.sdf
            WriteFragmentPDB(sidechain, '%s_frag%02d.pdb'%(basename.split(".")[0],count))
    #elif cut_type == "all":
    #   newrings, biggest_index_list = get_large_ring(cut_type, mol=mol,confId=-1)

    #   for biggest_index in biggest_index_list:
    #   # loop for different cores
    #       coreIdxs, sidechain_connectIdxs, core_connectIdxs = GetMacroCyclicFragmentIdx(newrings, biggest_index, mol=mol, confId=-1)
    #       core = GetFragmentConf(coreIdxs, core_connectIdxs, mol, confId=-1)
    #       count = 1
    #       WriteFragmentSDF(core, '%s_core%02d_frag%02d.sdf'%(basename.split(".")[0], biggest_index, count)) ## write down core.sdf
    #   
    #        for m, n in zip(sidechain_connectIdxs, core_connectIdxs):
    #            sidechainIdx = find_atom_ids_to_delete(mol, m, n)
    #            sidechain = GetFragmentConf(sidechainIdx, [m], mol)
    #            count = count + 1
    #            WriteFragmentSDF(sidechain, '%s_core%02d_frag%02d.sdf'%(basename.split(".")[0],biggest_index,count)) ## write down sidechain.sdf
#
    return None

def main(): 
    args = sys.argv[1:]

    if not args:
        print ('usage: python CyclicPoseFragmentation_noH.py [--input] filename.format')

        sys.exit(1)

    elif sys.argv[1] == '--help':
        print ('usage: python CyclicPoseFragmentation_noH.py [--input] filename.format')

        sys.exit(1)

    elif sys.argv[1] == '--input':
        infile = sys.argv[2]
        basename = os.path.basename(infile)
        print(basename)
        cut_type = "only1"
        if infile.split(".")[1] == "sdf":
            mol = Chem.SDMolSupplier(infile, removeHs=False)[0]
        elif infile.split(".")[1] == "mol2":
            print("mol2")
            mol = Chem.MolFromMol2File(infile, removeHs=False)
        runMethod(mol, basename, cut_type)
    

    else:
        sys.exit(1)

if __name__ == '__main__':
    main()




    

