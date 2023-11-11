'''
JWO 11/11/2023
'''
import sys
import numpy as np
from itertools import combinations 

num_to_char = {0:'A', 1:'C', 2:'G', 3:'T'}

# decimal -> kmer (base 4)
def dec_to_kmer(val, l):
    kmer = ['' for i in range(0, l)]
    for i in range(0, l):
        kmer[l - 1 - i] = num_to_char[val % 4]
        val //= 4
    return ''.join(kmer)


def main(argv = sys.argv):
    if(len(argv) != 4):
        print("usage: {0} {l (e.g. 11)} {k (e.g. 7)} {ofile}")
        sys.exit()

    l, k = int(argv[1]), int(argv[2])
    ofile = open(argv[3], 'w')


    combs = list(combinations(range(0, l), k)) # gapped-position combinations
    gkms = []
    for val in range(0, 4**k): # loop through all kmers
        kmer = dec_to_kmer(val, k)
        for comb in combs:
            gkm = ['-' for a in range(0, l)]
            for i, c in enumerate(comb):
                gkm[c] = kmer[i]
            gkms.append(''.join(gkm))
                

    for gkm in gkms:
        ofile.write(gkm + "\n")
    ofile.close()
main()
