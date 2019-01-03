import os

def get_num_fragments(infile):
    num = len([line for line in open(infile)])

    return num


def main(fn,outfile):
    #if fn in os.listdir("/Users/jianinglu1/Documents/PDBbind_v2016_refined/"):
        #datadir = "/Users/jianinglu1/Documents/PDBbind_v2016_refined/" + fn + "/SF_fragments_noH"
    #else:
        #datadir = "/Users/jianinglu1/Documents/general-set-except-refined/" + fn + "/SF_fragments_noH"
    #datadir = "/Users/jianinglu1/Documents/Train_Data/CSAR_decoy/CSAR-decoy/decoys-split/" + fn  + "_new/SF_fragments_noH"
    datadir = "/scratch/jl7003/BACE/prepare_final/pose_final/" + fn
    
    # CAS
    #datadir  = fn 

    #infile = datadir + "/Vina_frag.csv"
    infile = datadir + "/Vina_score_0.csv"
    num = get_num_fragments(infile)
    outfile.write(fn.split("/")[-1] + "," + str(num) + "\n")
    # CAS
    #outfile.write(fn.split("/")[-2] + "," + str(num) + "\n")

if __name__ == "__main__":
    outfile = open("../prepare_final/BACE_final_num_fragments.csv","w")
    outfile.write("pdb,num_frag\n")
    pdblist = ['{:0>3}'.format(i) for i in range(1,155)]
    for fn in pdblist:
        print(fn)
        main(fn,outfile)
    outfile.close()
