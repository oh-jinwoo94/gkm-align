'''
Jin Woo Oh

*Purpose: 
  LASTZ output
  "Negative strand intervals are measured from the 5' end of the query’s negative strand"
  For (-)strand matched query coordinates, flip coordinate to begin from coordinate=1.


* input file example 
  head -3 short_sequence_matches.axt
  mm10    chr17   66149955        66150097        hg38    chr1    12600   12742   +
  mm10    chr17   66150342        66150419        hg38    chr1    12985   13062   +

* specify query genome (e.g., hg38), use chrom.size to reverse coordinates for last column = "-". 
'''

import sys

def main(argv = sys.argv):
    if(len(argv) != 4):
        print("{0} {anchor file} {query genome (e.g. hg38)} {query genome chrom.sizes file}")
        sys.exit()

    qgenome = argv[2]
    chr_to_size = dict()
    with open(argv[3], 'r') as ifile:
        for line in ifile:
            words = line.split('\t')
            chrom, size = words[0], int(words[1])
            chr_to_size[chrom] = size 


    with open(argv[1]) as ifile:
        for line in ifile:
            words = line.rstrip().split()
            g_1, chrom_1, loc_1, g_2, chrom_2, loc_2, ori = words[0], words[1], (eval(words[2]) + eval(words[3])) // 2, words[4], words[5], (eval(words[6]) + eval(words[7])) // 2, words[8]

            if(ori == "-"): # lastz anchors has coordinate from the end of the genome for opposing strand
                if(g_2 == qgenome):
                    loc_2 = chr_to_size[chrom_2] - loc_2 - 1
                elif(g_1 == qgenome):
                    loc_1 = chr_to_size[chrom_1] - loc_1 - 1
                else:
                    print("query genome must be one of: " + g_1 + ", " + g_2)
                    sys.exit()
            print("\t".join([g_1,  chrom_1, str(loc_1), g_2, chrom_2, str(loc_2), ori]))
main()
