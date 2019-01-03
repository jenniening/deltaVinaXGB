import numpy as np
import get_pdbinfo




def get_RW(fn, inpro):
    '''
    Select the HOH in [2.0, 3.5] of protein

    '''

    outfile = open("test/RW_info.txt","w")

    pro = get_pdbinfo.pdbinfo(fn,file = inpro)
    pro_atoms = pro.getPolarAtoms()
    protein,waters = get_pdbinfo.pdbinfo(fn, lines = pro_atoms).getProteinWaters()
    waters_coord = get_pdbinfo.pdbinfo(fn,lines = waters).getCoords()
    protein_coord = get_pdbinfo.pdbinfo(fn,lines = protein).getCoords()
    ### calculate distance ###
    waters_coord = np.expand_dims(waters_coord, 1)
    protein_coord = np.expand_dims(protein_coord,0)
    distance = np.linalg.norm(waters_coord - protein_coord, axis = 2)
    distance_min = np.min(distance, axis = 1)
    rw_index = []
    for idx, i in enumerate(distance_min):
        if i >= 2.0 and i <= 3.5:
            rw_index.append(idx)
    for i in rw_index:
        rw_line = waters[i]
        rw_distance = distance[i]
        rw_chain = get_pdbinfo.chid(rw_line)
        if rw_chain != " ":
            rw_name = str(int(get_pdbinfo.resi(rw_line))) + "." + get_pdbinfo.chid(rw_line)
        else:
            rw_name = str(int(get_pdbinfo.resi(rw_line)))
            
        for idx, d in enumerate(rw_distance):
            if d <= 3.5 and d >= 2.0:
                pro_idx = idx
                pro_line = protein[idx]
                pro_chain = get_pdbinfo.chid(pro_line)
                pro_name = get_pdbinfo.resn(pro_line)
                if pro_chain != " ":
                    pro_idx = str(int(get_pdbinfo.resi(pro_line))) + "." + pro_chain 
                else:
                    pro_idx = str(int(get_pdbinfo.resi(pro_line)))
                pro_aname = get_pdbinfo.atmn(pro_line).strip()
                outline = fn + "," + pro_name + "," + pro_idx + "," + pro_aname + "," + rw_name + "," + str(round(d,2)) + "\n"
                outfile.write(outline)
    outfile.close()
    return None

if __name__ == "__main__":
    fn = "3c2f"
    inpro = "test/3c2f_protein_all_SF_type2.pdb"
    get_RW(fn, inpro)

