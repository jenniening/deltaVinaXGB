import os
import fileinput
import sys
import pandas as pd
import numpy as np
import get_pdbinfo
import get_inputtype
from get_inputtype import get_inputtype

if sys.platform == "linux":
    import software_path_linux as path
elif sys.platform == "darwin":
    import software_path_mac as path



MGLPY = path.path_mgl_python()
MGLUTIL = path.path_mgl_script()
vina = path.path_vina()

def get_coordinates_from_mol2(infile):
    '''
    Get coordinates from mol2 file

    '''
    init_lines = [line for line in open(infile)]
    atom_index = [index for index, line in enumerate(init_lines) if line[0:13] == "@<TRIPOS>ATOM"][0]
    bond_index = [index for index, line in enumerate(init_lines) if line[0:13] == "@<TRIPOS>BOND"][0]
    lines = [line for index, line in enumerate(init_lines)][atom_index + 1:bond_index]
    #coordinates = [ "%8s%8s%8s"%(line[16:46].split()[0][0:-1], line[16:46].split()[1][0:-1], line[16:46].split()[2][0:-1]) for line in lines]
    coordinates = [str("%8.3f%8.3f%8.3f"%(float(line[16:46].split()[0]), float(line[16:46].split()[1]),float(line[16:46].split()[2]))) for line in lines]

    return coordinates

def get_coordinates_from_pdb(infile):
    '''
    Get coordinates from pdb file

    '''

    lines = get_pdbinfo.pdbinfo(file = infile).getAtoms()
    Coords = get_pdbinfo.pdbinfo(lines = lines).getCoords()
    coordinates = [str("%8.3f%8.3f%8.3f"%(coord[0], coord[1],coord[2])) for coord in Coords]
    
    return coordinates


def get_ref_infor(inlig, ref_list, noname = True):
    '''
    Get reference coordinates and charge

    '''
    ref_infor = {"ref_" + str(num): [] for num in range(len([i for i in ref_list])) }
    if noname:
        # use lines order to finder order
        ref_lines = {"ref_" + str(num):None for num in range(len([i for i in ref_list])) }
        # here, we need the initial crystal mol2 file as the reference 
        ref_old = inlig
        if get_inputtype(ref_old) == "mol2":
            ref_old_lines = get_coordinates_from_mol2(ref_old)
        elif get_inputtype(ref_old) == "pdb":
            ref_old_lines = get_coordinates_from_pdb(ref_old)
        # get the lines for each file
        for num in range(len([i for i in ref_list])):
            if num == 0:
                ref_lines["ref_" + str(num)] = [line for line in open(ref_list[0]) if (line[0:4] == "ATOM") or (line[0:6] == "HETATM")]
            else:
                ref_lines["ref_" + str(num)] = [line[30:54] for line in open(ref_list[num]) if (line[0:4] == "ATOM") or (line[0:6] == "HETATM")]
        # get the new lines for ref_infor based on order from ref_old
        for index, line in enumerate(ref_old_lines): 
            write = False                                  
            for num in range(len([i for i in ref_list])):
                if num == 0:
                    if len([l for l in ref_lines["ref_0"] if line in l ]) != 0:
                        ref_line = [l for l in ref_lines["ref_0"] if line in l ][0]
                        write = True
                        ref_infor["ref_" + str(num)].append(ref_line)
                else:
                    if write:
                        new_line = ref_line[0:30] + ref_lines["ref_" + str(num)][index] + ref_line[54:]
                        ref_infor["ref_" + str(num)].append(new_line)
    
    else:
        # use atom name to find order, only for structure from PDBbind or CSAR
        for num in range(len([i for i in ref_list])):
            if num == 0:
                ref_infor["ref_" + str(num) ] = [line for line in open(ref_list[0]) if (line[0:4] == "ATOM") or (line[0:6] == "HETATM")]
            else:
                ref_infor["ref_" + str(num)] = []
                for line in ref_infor["ref_0"]:
                    name = line.split()[2]
                    ref_lines = [line for line in open(ref_list[num]) if (line[0:6] == "HETATM") or (line[0:4] == "ATOM")]
                    ref_infor["ref_" + str(num)].append([line for line in ref_lines if (line.split()[2] == name)][0])
    return ref_infor



def get_atom_index(ref_infor, frag_file, one_atom, atom_type):
    '''
    Get updated coordinates and charge for fragments
    
    '''
    new_lines = {"frag_" + str(num): [] for num in range(len(ref_infor.keys()))}
    lines_frag = [line for line in open(frag_file) if line[0:6] == "HETATM"]
    heavy_atom_length = len([line for line in lines_frag if line.split()[2][0:1] != "H"])
    coordinates_frag = [line[30:54] for line in lines_frag]
    n = 0 
    for idx, c_f in enumerate(coordinates_frag):

        # determine whether the frag is an atom 
        if one_atom and atom_type == "N3":
            for num in range(len(ref_infor.keys())):
                if num == 0:
                    [idx_prev, ref_line] = [[idx_r, line] for idx_r, line in enumerate(ref_infor["ref_" + str(num)]) if "%8s%8s%8s"%(c_f.split()[0],c_f.split()[1],c_f.split()[2]) in line][0]
                    new_lines["frag_" + str(num)].append("HETATM    1" + ref_line[11:])
                else:
                    new_lines["frag_" + str(num)].append("HETATM    1" + ref_infor["ref_" + str(num)][idx_prev][11:])

        elif one_atom and atom_type in ["One","CH"]:
            if lines_frag[idx].split()[2][0:1] != "H":
                n += 1
                if n <= heavy_atom_length:
                    for num in range(len(ref_infor.keys())):
                        if num == 0:
                            [idx_prev, ref_line] = [[idx_r, line] for idx_r, line in enumerate(ref_infor["ref_" + str(num)]) if "%8s%8s%8s"%(c_f.split()[0],c_f.split()[1],c_f.split()[2]) in line][0]
                            new_lines["frag_" + str(num)].append("HETATM%5s"%str(n) + ref_line[11:])
                        else:
                            new_lines["frag_" + str(num)].append("HETATM%5s"%str(n) + ref_infor["ref_" + str(num)][idx_prev][11:])
                else:
                    break

        else:
            for num in range(len(ref_infor.keys())):
                
                if num == 0:
                    [idx_prev, ref_line] = [[idx_r, line] for idx_r, line in enumerate(ref_infor["ref_" + str(num)]) if c_f in line][0]
                    new_lines["frag_" + str(num)].append("HETATM" + "%5s"%(str(idx + 1))+ ref_line[11:])
                else:
        
                    new_lines["frag_" + str(num)].append("HETATM" + "%5s"%(str(idx + 1)) + ref_infor["ref_" + str(num)][idx_prev][11:])
           
    return new_lines


 
def write_frag_pdbqt(prev_frag_file, new_frag, new_frag_file_list):
    '''
    Write new pdbqt file with corrected charge for fragments (C, Co, Crwo)

    '''
    keys = ["frag_" + str(num) for num in range(len(new_frag_file_list))]

    if prev_frag_file != None:
        for idx, filename in enumerate(new_frag_file_list):
            n = 0 
            outfile = open(filename, "w")
            for line in open(prev_frag_file):
                if line[0:6] != "HETATM":
                    outfile.write(line)
                else:
                    outfile.write(new_frag[keys[idx]][n])
                    n +=1
            outfile.close()
    else:
        for idx, filename in enumerate(new_frag_file_list):
            outfile = open(filename, "w")
            outfile.write("REMARK  0 active torsions:\n"
                        + "REMARK  status: ('A' for Active; 'I' for Inactive)\n"
                        + "ROOT\n")
            for line in new_frag[keys[idx]]:
                outfile.write(line)
            outfile.write("ENDROOT\n"
                        + "TORSDOF 0\n")
            outfile.close()

    return None

     
def get_pdbqt(frag,prev_frag_pdbqt, datadir):
    '''
    Generate pdbqt file for fragments

    '''

    cmd = MGLPY + " "  + MGLUTIL + "prepare_ligand4.py -l " + frag + " -o " + prev_frag_pdbqt +  " -U 'nphs'"
    os.system(cmd)

    return None

def preparelig(inlig, ligpdbqt):
    """
    Prepare ligand PDBQT file by MGLTools 
    
    """
    print("Generate ligand pdbqt")

    cmd = MGLPY + " "  + MGLUTIL + "prepare_ligand4.py -l " + inlig  + " -o " + ligpdbqt +  " -U 'nphs'"
    os.system(cmd)

def prepareprot(inprot, protpdbqt):
    """
    Prepare ligand PDBQT file by MGLTools 
    
    """
    print("Generate protein pdbqt")

    cmd = MGLPY + " "  + MGLUTIL + "prepare_receptor4.py -r " + inprot  + " -o " + protpdbqt +  " -U 'nphs'" 
    os.system(cmd)


def get_ref(fn, inlig, datadir, min = False, min_RW = False, RW = False, min_BW = False, min_PW = False, decoy = False, decoy_list = None):
    '''
    Get reference file

    '''
    ref_list = []
    if not os.path.isfile(os.path.join(datadir,fn + "_ligand_Vina58.pdbqt")):
        inlig_out = os.path.join(datadir,fn +"_ligand_Vina58.pdbqt")
        preparelig(inlig, inlig_out)

    ref_list.append(os.path.join(datadir,fn +"_ligand_Vina58.pdbqt"))
    if min:
        ref_list.append(os.path.join(datadir,fn + "_lig_min.pdb"))
    if min_RW:
        ref_list.append(os.path.join(datadir,fn +"_lig_min_RW.pdb"))
    if RW:
        ref_list.append(os.path.join(datadir,fn +"_ligand_Vina58.pdbqt"))
    if min_BW:
        ref_list.append(os.path.join(datadir,fn + "_lig_min_BW.pdb"))
    if min_PW:
        ref_list.append(os.path.join(datadir,fn + "_lig_min_PW.pdb"))
    if decoy:
        for decoy in decoy_list:
            ref_list.append(os.path.join(datadir, decoy))
    return ref_list

def get_prot(fn,datadir, min = False, min_RW = False, RW = False, min_BW = False, min_PW = False,decoy = False, decoy_list = None, decoy_pro = None):
    '''
    Get protein pdbqt

    '''
    prot_list = []
    if not os.path.isfile(os.path.join(datadir, fn +"_protein_Vina58.pdbqt")):
        inprot = os.path.join(datadir,fn +"_protein.pdb")
        inprot_out = os.path.join(datadir,fn + "_protein_Vina58.pdbqt")
        prepareprot(inprot,inprot_out)
    prot_list.append(os.path.join(datadir, fn + "_protein_Vina58.pdbqt"))
    if min:
        prot_list.append(os.path.join(datadir, fn + "_protein_Vina58.pdbqt"))
    if min_RW:
        if not os.path.isfile(os.path.join(datadir, fn +"_protein_RW_Vina58.pdbqt")):
            inprot = os.path.join(datadir,fn +"_protein_RW.pdb")
            inprot_out = os.path.join(datadir,fn +"_protein_RW_Vina58.pdbqt")
            prepareprot(inprot,inprot_out)
        prot_list.append(os.path.join(datadir, fn + "_protein_RW_Vina58.pdbqt"))
    if RW:
        prot_list.append(os.path.join(datadir, fn + "_protein_RW_Vina58.pdbqt"))
    if min_BW:
        if not os.path.isfile(os.path.join(datadir, fn +"_protein_BW_Vina58.pdbqt")):
            inprot = os.path.join(datadir,fn +"_protein_BW.pdb")
            inprot_out = os.path.join(datadir,fn +"_protein_BW_Vina58.pdbqt")
            prepareprot(inprot,inprot_out)
        prot_list.append(os.path.join(datadir, fn + "_protein_BW_Vina58.pdbqt"))
    if min_PW:
        if not os.path.isfile(os.path.join(datadir, fn +"_protein_PW_Vina58.pdbqt")):
            inprot = os.path.join(datadir,fn +"_protein_PW.pdb")
            inprot_out = os.path.join(datadir,fn +"_protein_PW_Vina58.pdbqt")
            prepareprot(inprot,inprot_out)
        prot_list.append(os.path.join(datadir, fn + "_protein_PW_Vina58.pdbqt"))

    if decoy:
        if not os.path.isfile(os.path.join(datadir,decoy_pro.split(".")[0] + "_Vina58.pdbqt")):
            inprot = os.path.join(datadir,decoy_pro)
            inprot_out = os.path.join(datadir,decoy_pro.split(".")[0] + "_Vina58.pdbqt")
            prepareprot(inprot,inprot_out)
        num = len(decoy_list)
        for i in range(num):
            prot_list.append(os.path.join(datadir,decoy_pro.split(".")[0] + "_Vina58.pdbqt"))


    return prot_list


def get_frag(fn,frag_dir):
    frag_num =  len([ filename for filename in os.listdir(frag_dir) if ("frag" in filename) and (filename.endswith(".pdb"))])
    frag_list = [os.path.join(frag_dir,'%s_ligand_frag%02d.pdb'%(fn, i) )for i in range(1,frag_num + 1)]


    return frag_list



def runVina(fn, frag_id, protpdbqt, ligpdbqt, datadir ):
    """Run modified AutoDock Vina program with Vina score and 58 features """

    frag_id +=1
    cmd = vina + " --receptor " + protpdbqt + " --ligand " + ligpdbqt + "  --score_only --log " + os.path.join(datadir,fn + "_" + str(frag_id) + "_score.txt") + " > " + datadir + "/out_vina.log"
    os.system(cmd)

    vinalist = [fn,str(frag_id)]
    for lines in fileinput.input(os.path.join(datadir, fn + "_" + str(frag_id) + "_score.txt")):
        if lines[0:4] in ["Affi", "Term"]:
            if sys.platform == "linux":
                ### vina in linux ###
                line = lines.strip("\n").split()
            elif sys.platform == "darwin":
                ### vina in os ###
                line = lines.strip("\n").split(":")
            vinalist.append(line[1])

    if len(vinalist) != 61:
        vinalist = [fn,str(frag_id)] + ['NA' for i in range(59)]

    return vinalist


    

def run_Vina_Fragment(fn, inlig, datadir, datadir_frag, min = False, min_RW = False, RW = False, min_BW = False, min_PW = False, decoy = False, decoy_list = None, decoy_pro = None):
    
    inlig = os.path.join(datadir, inlig)
    # get ref_lig name
    ref_list = get_ref(fn,inlig, datadir,min, min_RW, RW, min_BW, min_PW, decoy, decoy_list)
    # get protein pdbqt
    prot_list = get_prot(fn,datadir,min,min_RW, RW, min_BW, min_PW, decoy, decoy_list, decoy_pro)
    # get previous ligand pdbqt 
    ref_infor = get_ref_infor(inlig, ref_list,noname = True)
    # get fragments name
    frag_dir = datadir_frag
    frag_list = get_frag(fn, frag_dir)
    # newfile names
    outfile_list = [open(frag_dir + "/Vina_score_" + str(num) + ".csv","w") for num in range(len(ref_list))]

    for frag_id, frag in enumerate(frag_list):
        # if fragment has no polar atom or just one atom, it can't generate pdbqt file without error message
        one_atom = False
        atom_type = None
        lines = [line for line in open(frag) if line[0:6] == "HETATM"]
        if len(lines) == 1:
            one_atom = True
            atom_type = "One"
        elif len(lines) == 4 and len([line for line in lines if line.split()[2][0:1] == "H"]) == 3 and len([line for line in lines if line.split()[2][0:1] == "C"]) == 1:
            ### CH3 ###
            one_atom = True
            atom_type = "CH"
        elif len(lines) == 6 and len([line for line in lines if line.split()[2][0:1] == "H"]) == 3 and len([line for line in lines if line.split()[2][0:1] == "C"]) == 3:
            ### CCCH3 ###
            one_atom = True
            atom_type = "CH"
        elif len(lines) == 3 and len([line for line in lines if line.split()[2][0:1] == "N"]) == 3:
            ### NNN ###
            one_atom = True
            atom_type = "N3"
        elif len(lines) == 3 and len([line for line in lines if line.split()[2][0:1] == "N"]) == 1 and len([line for line in lines if line.split()[2][0:1] == "C"]) == 1 and len([line for line in lines if line.split()[2][0:1] == "S"]) == 1 and fn == "3ary":
            ### SCN in 3ary ###
            one_atom = True
            atom_type = "N3"
        if one_atom:
            prev_frag_file = None
            new_lines= get_atom_index(ref_infor, frag, one_atom, atom_type)
        else:
            prev_frag_file = frag.split(".")[0] + ".pdbqt"
            if os.path.isfile(prev_frag_file) == False:
                # generate pdbqt file for fragment
                get_pdbqt(frag, prev_frag_file, frag_dir)
            new_lines= get_atom_index(ref_infor, prev_frag_file, one_atom, atom_type)
        # correct charge and generate new pdbqt file 
        total_files = len([i for i in ref_list])
        new_frag_list = [frag.split(".")[0] + "_correct_" + str(num) + ".pdbqt" for num in range(total_files)]
        write_frag_pdbqt(prev_frag_file, new_lines, new_frag_list)
        # calculate Vina58 features
        for idx, new_frag in enumerate(new_frag_list):
            vina_list = runVina(fn, frag_id, prot_list[idx], new_frag, frag_dir)
            outfile_list[idx].write(",".join(vina_list) + "\n")
    for outfile in outfile_list:
        outfile.close()
    return None


def generate_data(fn, data_type,datadir):
    '''
    Generate core/side Vina Score 
    '''
    ### data_type = 0 --> cry; data_type = 1 ---> Co; data_type = 2 --> Crwo ###
    columns = ["pdb","idx","vina"] + ["vina" + str(i) for i in range(1,59)]
    infile = os.path.join(datadir, "Vina_score_" + data_type + ".csv")
    if len([line for line in open(infile)]) == 1:
        ### no sidechain and whole structure is core ###
        num_frag = 0
    else:
        num_frag = len([line for line in open(infile)])

    if os.stat(infile).st_size != 0:

        df = pd.read_csv(infile, header = None, dtype = {0:str})
        df.columns = columns

        df["vina"] = df["vina"] * -0.73349
        df["pdb"] = df["pdb"].astype(str)

        ### core Vina score is the first one
        df_core = df[df["idx"] == 1].copy()
        df_core.drop("idx",axis=1,inplace=True)
        df_core.to_csv(os.path.join(datadir,"Vina_core_" + data_type + ".csv"), index = False)
        if num_frag != 0:
            ### sum scores for all side fragments ###
            df_side = df[df["idx"] != 1].copy().groupby("pdb").agg({col:np.sum for col in df.columns[2:]})
            df_side.reset_index(inplace = True)
            df_side.to_csv(os.path.join(datadir,"Vina_side_" + data_type + ".csv"), index = False)
        else:
            df_side = open(os.path.join(datadir,"Vina_side_" + data_type + ".csv"),"w")
            df_side.write(",".join(["pdb","vina"] + ["vina" + str(i) for i in range(1,59)]) + "\n")
            df_side.write(",".join([fn] + ["0" for i in range(59)]) + "\n")
            df_side.close()
    else:
        sys.exit("Error: No Vina Score for Fragments")

    return num_frag



def combine_data(fn,datadir,data_type, outfile_core, outfile_side):
    '''
    combine data for pdblist
    '''


    num_frag = generate_data(fn,data_type,datadir)
    ### write out for many structures ###
    for line in open(os.path.join(datadir,"Vina_core_" + data_type + ".csv")):
        if line.split(",")[0] != "pdb":
            # For CatS--> frag name is not same as pdb name
            outfile_core.write(fn + "," + ",".join(line.split(",")[1:]))
            #outfile_core.write(line)
    for line in open(os.path.join(datadir, "Vina_side_" + data_type + ".csv")):
        if line.split(",")[0] != "pdb":
            # For CatS--> frag name is not same as pdb name
            outfile_side.write(fn + "," + ",".join(line.split(",")[1:]))
            #outfile_side.write(line)        

    return num_frag
    
def main():
    args = sys.argv[1:]

    if not args:
        print ('usage: python Get_Fragment_Vina58_specific.py [--pdbid] fn [--datadir] dir [--fragdir] dir')
        sys.exit(1)
    elif sys.argv[1] == '--help':
        print ('usage: python Get_Fragment_Vina58_specific.py [--pdbid] fn [--datadir] dir [--fragdir] dir')
        sys.exit(1)

    elif sys.argv[1] == '--pdbid':
        fn = sys.argv[2]
        if sys.argv[3] == "--datadir":
            datadir = sys.argv[4]
            if sys.argv[5] == "--fragdir":
                fragdir = sys.argv[6]
            else:
                fragdir = sys.argv[4]
                #run_Vina_Fragment(fn,datadir,fragdir)
    else:
        sys.exit(1)


if __name__ == "__main__":
    #### ------------BACE_1 ----------------------------
    datadir  = "/scratch/jl7003/BACE/prepare_final/"
    fragdir = "/scratch/jl7003/BACE/prepare_final/pose_final/"
    pdblist = ['{:0>3}'.format(i) for i in range(1,155)]
    data_type = "3"
    outfile_core = open(datadir + "BACE_final_core_" + data_type + "_new.csv","w")
    outfile_core.write("pdb,vina," + ",".join(["vina" + str(i) for i in range(1,59)]) + "\n")
    outfile_side = open(datadir + "BACE_final_side_" + data_type + "_new.csv","w")
    outfile_side.write("pdb,vina," + ",".join(["vina" + str(i) for i in range(1,59)]) + "\n")
    ##previous_file = "/Users/jianinglu1/Documents/prepare_01/Users/jianinglu1/Documents/prepare03/Input_min_RW_prepare_03.csv"
    for fn in pdblist:
        print(fn)
        datadir_current  = datadir + fn  + "_new"
        os.system("cp -r " + fragdir + fn + " " + fragdir + fn + "_new")
        fragdir_current = fragdir + fn + "_new"
        #run_Vina_Fragment(fn,datadir_current,fragdir_current, True, True, True, False)

        combine_data(fn,fragdir_current, data_type, outfile_core, outfile_side)
    outfile_core.close()
    outfile_side.close()
    #### -------------CAS -------------------------
   #datadir  = "/Users/jianinglu1/Documents/ligand_align_core_gc3/"
   #datadir_list = [ datadir + file for file in os.listdir(datadir) if "CatS" in file]
   #data_type = "1"
   #previous_file = None
   #outfile_core = open(datadir + "CAS_core_core_" + data_type + ".csv","w")
   #outfile_core.write("pdb,vina," + ",".join(["vina" + str(i) for i in range(1,59)]) + "\n")
   #outfile_side = open(datadir + "CAS_core_side_" + data_type + ".csv","w") 
   #outfile_side.write("pdb,vina," + ",".join(["vina" + str(i) for i in range(1,59)]) + "\n")
   #for idx, datadir in enumerate(datadir_list):
   #    print(datadir)
   #    fragdir_current = datadir + "/fragments"
   #    fn = [filename for filename in os.listdir(datadir) if filename.startswith("conformer")][0]
   #    print(idx, fn)
   #    #datadir = datadir + "/" + fn 
   #    #run_Vina_Fragment(fn,datadir,fragdir_current, True, False)
   #    combine_data(datadir.split("/")[-1],fragdir_current, data_type, previous_file, outfile_core, outfile_side)
   #outfile_core.close()
   #outfile_side.close()
   #    



