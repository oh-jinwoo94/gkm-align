import pandas as pd
import sys

def main(argv= sys.argv):
    if(len(argv) != 3):
        print("usage: <0> <list of weight files to avgerage> <ofile>")
        sys.exit(1)

    lfile = open(argv[1], 'r')
    ofile = open(argv[2], 'w')
    n = 0
    for i, line in enumerate(lfile):
        n += 1

        fname = line.rstrip()
        table = pd.read_csv(fname, sep='\t', header=None) 
        if(i==0):
            print(table.values)
            kmers = table.values[:,0]
            w = table.values[:,1]
            print(w[0])
        else:
            w += table.values[:,1]
            print(w[0])
    if(n==0):
        print("no file in the list file")
        sys.exit(1)
    else:
        w /= n
        print(w[0])



    for i, kmer in enumerate(kmers):
        ofile.write("\t".join([kmer, str(w[i])]) + '\n')

    ofile.close()
    lfile.close()   

main()
