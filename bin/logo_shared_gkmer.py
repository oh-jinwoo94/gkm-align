# For generating gkm-align mapping visuaslization (e.g, Fig 5C of Oh & Beer Nature comm 2024)
# Upon request on "Issues", I will give more details on how to run this script. 

import logomaker
import pandas as pd
import sys
from matplotlib import pyplot as plt

def main(argv = sys.argv):
    
    if(len(argv)!=7):
        print("Usage <0> <matrix file 1> <matrix file 2> <ID for matrix file 1> <ID for matrix file 2> <info> <ofile name>")
        sys.exit(0)
    fig, (ax_q, ax_t) = plt.subplots(nrows=2)

    mfile1 = argv[1]
    mfile2 = argv[2]
    id1 = argv[3]
    id2 = argv[4]
    cons_info = argv[5]
    
    ofname = argv[-1]
    mat_q = pd.read_csv(mfile1, delim_whitespace=True, index_col=0, comment='#')
    logo_q = logomaker.Logo(mat_q, ax = ax_q, center_values=False, figsize = [100, 5])

    mat_t = pd.read_csv(mfile2, delim_whitespace=True, index_col=0, comment='#')
    logo_t = logomaker.Logo(mat_t, ax = ax_t, center_values=False, figsize = [100, 5])

    ax_q.text(0.1, 0.9, cons_info)
    ax_q.set_title(id1)
    ax_t.set_title(id2)

    ax_q.set_xlabel("position")
    ax_t.set_xlabel("position")

    ax_q.set_ylabel("gkm cons.")
    ax_t.set_ylabel("gkm cons.")

    fig.set_size_inches(30, 4)
    fig.tight_layout()
    fig.savefig(ofname)

main()
