'''
simple smith-waterman for mapping interpretation 

.mapped file format 
with no header
hg38    chr18   22084683        22085138        -->     mm10    chr18   10986037        10986491
'''
import numpy as np
import sys
np.set_printoptions(threshold = np.inf)
sys.setrecursionlimit(5000)

def edit_distance(align1, align2):
    indel = len([a for a in align1 if a=='-'] + \
                [a for a in align2 if a=='-'])
    sub = len([i for i in range(0, len(align1)) if \
                (align1[i]!='-' and align2[i]!='-') and \
                (align1[i]!=align2[i])])
    return(indel+sub)

def build_D(alpha, beta, d_sub, d_indel):
    # initialize dynamic matrix
    D = np.zeros([len(alpha)+1, len(beta)+1])
    for i in range(1,D.shape[0]):
        D[i,0] = i*d_indel
    for j in range(1,D.shape[1]):
        D[0,j] = j*d_indel

    # fill in dynamic matrix
    for i in range(1, D.shape[0]):
        for j in range(1, D.shape[1]):
            if(alpha[i-1] == beta[j-1]):
                cand1 = D[i-1,j-1] 
            else:
                cand1 = D[i-1,j-1] + d_sub
            cand2 = D[i-1,j] + d_indel
            cand3 = D[i,j-1] + d_indel
            D[i,j] = min([cand1, cand2, cand3])
    return D

def D_deep_search(D, i, j, alpha, beta, align1, align2, d_sub, d_indel):
    if(not(i==0 and j==0)):
        #Boolean values that encode which 'directions' are possible
        match = ((i>0 and j>0) and  ((D[i-1,j-1])==D[i,j]) and \
                (alpha[i-1]==beta[j-1]))
        sub = ((i>0 and j>0) and ((D[i-1,j-1] + d_sub) == D[i,j]) and \
                (alpha[i-1]!=beta[j-1]))
        indel1 = ((i>0) and  ((D[i-1,j] + d_indel) == D[i,j]))
        indel2 = ((j>0) and ((D[i,j-1] + d_indel) == D[i,j]))

        if(match or sub):
            D_deep_search(D, i-1, j-1, alpha, beta,\
                    alpha[i-1] + align1, beta[j-1] + align2, d_sub, d_indel)
        if(indel1):
            D_deep_search(D, i-1, j, alpha, beta, \
                    alpha[i-1] + align1,'-' + align2, d_sub, d_indel)
        if(indel2):
            D_deep_search(D, i, j-1, alpha, beta, \
                    '-' + align1, beta[j-1] +align2, d_sub, d_indel)

    else: #maximum depth
        print(">" + align1)
        print(">" + align2)
        sys.exit()
def main(argv=sys.argv):
   
    if(len(argv) != 5):
        print("{0} {.mapped} {genome loc} {substitution penalty} {indel penalty}")
        sys.exit()
        
    ifile = open(argv[1], 'r')
    gdir = argv[2]
    d_sub, d_indel= eval(argv[3]), eval(argv[4])
    
  #hg38    chr18   22084683        22085138        -->     mm10    chr18   10986037        10986491   
    for line in ifile:
        print(line.rstrip())
        # read in the sequences
        words = line.split()
        fafile1 = open(gdir + "/" + words[0] + "/" + words[1] + ".fa", 'r')
        fafile1.readline()
        fa_seq1 = "".join([s.rstrip() for s in fafile1])
        begin1, end1 = int(words[2]), int(words[3])
        seq1 = fa_seq1[begin1:end1+1]

        fafile2 = open(gdir + "/" + words[5] + "/" + words[6] + ".fa", 'r')
        fafile2.readline()
        fa_seq2 = "".join([s.rstrip() for s in fafile2])
        begin2, end2 = int(words[7]), int(words[8])
        seq2 = fa_seq2[begin2:end2+1]


        seq1, seq2 = seq1.upper(), seq2.upper()
        D = build_D(seq1, seq2, d_sub, d_indel)
        (i,j) = D.shape
        D_deep_search(D, i-1, j-1, seq1, seq2, '','', d_sub, d_indel)
main()
