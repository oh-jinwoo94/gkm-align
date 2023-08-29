'''
filter anchors for those that lie within hg38-mm10 syntenic intergenic loci. 
(base) [joh27@troctolite syntenic_intergenic]$ head -3  /mnt/data0/joh27/projects/alignment_enhancer_conservation/analysis/gene_syntenny_and_line/hg38_mm10_ortholog_syntenic_intergenic_pairs.txt
species_1_genome_build  species_1_gene_name_1   species_1_gene_1_strand species_1_gene_name_2   species_1_gene_2_strand species_1_chr   species_1_range_begin    species_1_range_end     species_2_genome_build  species_2_gene_name_1   species_2_gene_1_strand species_2_gene_name_2   species_2_gene_2_strand species_2_chr    species_2_range_begin   species_2_range_end     m_genepair_orientation_relative_to_h
hg38    ENSG00000107099 +       ENSG00000107104 +       chr9    214854  746106  mm10    ENSMUSG00000052085      +       ENSMUSG00000032702      +       chr19    24999534        25434496        not-inverted
hg38    ENSG00000107104 +       ENSG00000137090 +       chr9    470291  969090  mm10    ENSMUSG00000032702      +       ENSMUSG00000024837      +       chr19    25236975        25604329        not-inverted

(base) [joh27@troctolite syntenic_intergenic]$ head anchors.axt
# target / query 
0 chr17 66149955 66150097 chr1 12600 12742 + 6237
1 chr17 66150342 66150419 chr1 12985 13062 + 3424
2 chr17 66150518 66150732 chr1 13165 13379 + 10391
3 chr2 111478783 111478958 chr1 52449 52624 + 6294
4 chr2 111292269 111292442 chr1 52451 52624 + 4174

'''

import sys

def main(argv = sys.argv):
    if(len(argv) != 5):
        print("{0} {anchor file} {syntenic intergenic file} {query genome (e.g. hg38)} {query genome chrom.sizes file}")
        sys.exit()

    qgenome = argv[-2]
    chr_to_size = dict()
    with open(argv[-1], 'r') as ifile:
        for line in ifile:
            words = line.split('\t')
            chrom, size = words[0], int(words[1])
            chr_to_size[chrom] = size 


    syn_interg_map = dict() # tchrom-qchrom -> tcoord-qcoord
    with open(argv[2]) as ifile:
        ifile.readline()
        for line in ifile:
            words = line.split('\t')
            if(words[0] == qgenome):
                chrom_key = words[13] + "-" + words[5]
                coord_range_pair = [int(words[14]), int(words[15]), int(words[6]), int(words[7])]
            else:
                chrom_key = words[5] + "-" + words[13]
                coord_range_pair = [int(words[6]), int(words[7]), int(words[14]), int(words[15])]

            try:
                syn_interg_map[chrom_key].append(coord_range_pair)
            except KeyError:
                syn_interg_map[chrom_key] = [coord_range_pair]
    with open(argv[1]) as ifile:
        for line in ifile:
            words = line.split()
            tchrom, tloc, qchrom, qloc, ori = words[1], (eval(words[2]) + eval(words[3])) / 2, words[4], (eval(words[5]) + eval(words[6])) / 2,\
                                              words[7]
            if(ori == "-"): # lastz anchors has coordinate from the end of the genome for opposing strand 
                qloc = chr_to_size[qchrom] - qloc - 1
                
            try:
                coord_range_pair_list = syn_interg_map[tchrom + "-" + qchrom]
                for crange_pair in coord_range_pair_list:
                    if(crange_pair[0] < tloc and tloc < crange_pair[1] and crange_pair[2] < qloc and qloc < crange_pair[3]):
                        print("\t".join(words))
                        break
            except KeyError:
                continue 
main()
