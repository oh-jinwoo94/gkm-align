'''
##
Input: 
1)
A single-fasta file 
>example
GATTAACTCCG

2) List of all gapped-kmers(l,k)
AAAAAAA----
AAAAAA-A---
AAAAAA--A--
.
.
##
Output:
AAAAAAA----     0
AAAAAA-A---     0
AAAAAA--A--     0
'''
import sys
import numpy as np
from itertools import combinations 

comp={'A':'T', 'T':'A', 'G':'C', 'C':'G'}
num_to_char = {0:'A', 1:'C', 2:'G', 3:'T', 4:'-'}

# returns reverse complment
def revcomp(seq):
    oseq=''
    for nc in seq:
        oseq=comp[nc]+oseq
    return np.array(list(oseq))

# kmer -> [gkm1, gkm2, ..]
def kmer_to_gkms(kmer, l, k):
    # choose (l-k) positions out of l position and replace with '-'
    combs = list(combinations(range(0, l), l-k))
    gkms = [] 
    for i, comb in enumerate(combs):
        gkmer = np.copy(np.array(list(kmer)))
        gkmer[list(comb)] = '-'
        gkms.append(''.join(gkmer))
    return gkms


def main(argv = sys.argv):
    if(len(argv) != 7):
        print("usage: {0} {single .fa} {gapped_kmer_list.txt} {gapped-kmer size (e.g. l=11)} {#non-gapped positions (e.g. k=7)} {include reverse complement? (y or n)} {ofile}")
        sys.exit()


    # read sequence
    ifile = open(argv[1], 'r')
    ifile.readline()
    seq = ifile.readline().rstrip()

    # gkm -> count
    gkm_to_cnt = {line.rstrip() : 0 for line in open(argv[2], 'r')}

     
    l, k = int(argv[3]), int(argv[4])
    rc = argv[5]
    ofile = open(argv[6], 'w')


    for i in range(0, len(seq)-l+1):
        kmer = seq[i:i+l]
        if(rc == 'y'):
            gkmers = kmer_to_gkms(kmer, l, k) + kmer_to_gkms(revcomp(kmer), l, k)
        elif(rc == 'n'):
            gkmers = kmer_to_gkms(kmer, l, k)
        else:
            print("y or n must be chosen")
            
        for gkmer in gkmers:
            gkm_to_cnt[gkmer] += 1

    for gkm in list(gkm_to_cnt):
        ofile.write("\t".join([gkm, str(gkm_to_cnt[gkm])]) + '\n')
    ofile.close()
main()
