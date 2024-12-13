#sample any genomic regions of given size

import os, sys
import random

k=0    #kmer table kmer length. Updated in main()
kmer_table=dict()

def main(argv=sys.argv):
    if len(argv)!=5:
        print("usage: {0}  <ofile> <length> <# seq> <random seed>")
        print(argv)
        sys.exit(0)
    print("Generating random sequences.")
    ofile = open(argv[1], 'w')
    length = int(argv[2])
    n = int(argv[3])
    rseed = int(argv[4])
    random.seed(rseed)

    for i in range(0,n):
        seq = "".join([random.choice(['A','T','G','C']) for j in range(0,length)])
        ofile.write('>rand_gen'+str(i)+'\n')
        ofile.write(seq+'\n')
    ofile.close()
     
    
main()

