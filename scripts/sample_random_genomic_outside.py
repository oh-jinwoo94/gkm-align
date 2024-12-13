# warning: uses >20 gigabites of RAM

# sample any genomic regions of given size
# takes extra bed input, so that sampling occurs outside the bed region
# intended for sampling outside of union(all DHS, TF peaks)

import os, sys
import random

k=0    #kmer table kmer length. Updated in main()
kmer_table=dict()

def main(argv=sys.argv):
    if len(argv)!=7:
        print("usage: {0} <genome directory> <.bed>  <ofile> <length> <# seq> <random seed>")
        print(argv)
        sys.exit(0)

    print("Sampling genomic sequences.")

    genome_dir=argv[1]
    bedfile_name = argv[2] 
    ofile_name = argv[3]
    length = int(argv[4])
    n = int(argv[5])
    rseed = int(argv[6])
    random.seed(rseed)


    # chr -> peak beds in the chromosome
    chrom_peaks = dict()
    if(bedfile_name not in ["none", "None", "NONE"]):
        with open(bedfile_name, 'r') as bedfile:
            for line in bedfile:
                words = line.split()
                chrom = words[0]; pos = [int(words[1]), int(words[2])]
                if(chrom not in chrom_peaks.keys()):
                    chrom_peaks[chrom] = [pos]
                else:
                    chrom_peaks[chrom].append(pos)
    else:
        pass # chrom_peaks will just remain an empty dictionary

    f_name_list = os.listdir(genome_dir)
    f_name_list = [name for name in f_name_list if ((name[-3:]=='.fa') and  ('_' not in name[name.index('chr')+3:name.index('.')]))]
    f_name_list = [name for name in f_name_list if name!="chrM.fa"]
    whole_genome = []
    for f_name in f_name_list:
        print(f_name)
        f = open(genome_dir + '/'+f_name, 'r')
        whole_chrom = []
        chrom = f_name[:-3]
        try:
            peaks = chrom_peaks[chrom]
        except KeyError:
            peaks = []
      
        for line in f: #read in the entire chromosome into a file 
            if(line[0]!='>'):
                whole_chrom += list(line.rstrip())

        for peak in peaks: #remove peak 
            for i in range(peak[0], peak[1]+1):
                whole_chrom[i] = '-'

        for i,nuc in enumerate(whole_chrom):
            if(nuc == 'n' or nuc == 'N'):
                whole_chrom[i] = '-'
        whole_genome += whole_chrom
        f.close()

    print("Finished reading genome")
    ofile = open(ofile_name,'w')

    n_sampled = 0
    while(n_sampled < n): #sample until target sample size reached
        if(n_sampled%100 == 0):
            print(str(n_sampled) + '/' + str(n))

        start = random.randint(0, len(whole_genome) - length)   
        if('-' not in whole_genome[start:start+length]):
            seq = whole_genome[start:start+length]
            ofile.write('>samp'+str(n_sampled)+'\n')
            ofile.write("".join(whole_genome[start:start+length]) +'\n')
            for i in range(start, start+length): # for preventing duplicate sequences
                whole_genome[i] = '-'
            n_sampled += 1 
    ofile.close()
main()
