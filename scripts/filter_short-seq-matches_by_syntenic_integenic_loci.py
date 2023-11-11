'''
Jin Woo Oh

two input files (1) short sequence match coordinates (2) human-mouse syntenic intergenic loci

(1)   
(base) [joh27@troctolite whole_genome]$ head -2 short_sequence_matches_coordinate_fixed.axt
mm10    chr17   66150026        hg38    chr1    12671   +
mm10    chr17   66150380        hg38    chr1    13023   +

// coordinates starting from 1 for both (+/-)


(2)
filter anchors for those that lie within hg38-mm10 syntenic intergenic loci. 
(base) [joh27@troctolite syntenic_intergenic]$ head -3  /mnt/data0/joh27/projects/alignment_enhancer_conservation/analysis/gene_syntenny_and_line/hg38_mm10_ortholog_syntenic_intergenic_pairs.txt
species_1_genome_build  species_1_gene_name_1   species_1_gene_1_strand species_1_gene_name_2   species_1_gene_2_strand species_1_chr   species_1_range_begin    species_1_range_end     species_2_genome_build  species_2_gene_name_1   species_2_gene_1_strand species_2_gene_name_2   species_2_gene_2_strand species_2_chr    species_2_range_begin   species_2_range_end     m_genepair_orientation_relative_to_h
hg38    ENSG00000107099 +       ENSG00000107104 +       chr9    214854  746106  mm10    ENSMUSG00000052085      +       ENSMUSG00000032702      +       chr19    24999534        25434496        not-inverted
hg38    ENSG00000107104 +       ENSG00000137090 +       chr9    470291  969090  mm10    ENSMUSG00000032702      +       ENSMUSG00000024837      +       chr19    25236975        25604329        not-inverted




'''

import sys

def main(argv = sys.argv):
    if(len(argv) != 3):
        print("{0} {short_seq_match file} {syntenic intergenic file}")
        sys.exit()


    syn_interg_map = dict()
    genome_1 = ""
    with open(argv[2]) as ifile:
        ifile.readline()
        for line in ifile:
            words = line.split('\t')
            genome_1 = words[0]
            chrom_key = words[5] + "-" + words[13]
            coord_range_pair = [int(words[6]), int(words[7]), int(words[14]), int(words[15])]

            try:
                syn_interg_map[chrom_key].append(coord_range_pair)
            except KeyError:
                syn_interg_map[chrom_key] = [coord_range_pair]
    with open(argv[1]) as ifile:
        for line in ifile:
            words = line.split()
            g_1, chrom_1, loc_1, g_2, chrom_2, loc_2 = words[0], words[1], eval(words[2]), words[3], words[4], eval(words[5])

            if(genome_1 == g_2):
                g_1, chrom_1, loc_1, g_2, chrom_2, loc_2 = g_2, chrom_2, loc_2, g_1, chrom_1, loc_1
                
            try:
                coord_range_pair_list = syn_interg_map[chrom_1 + "-" + chrom_2]
                for crange_pair in coord_range_pair_list:
                    if(crange_pair[0] < loc_1 and loc_1 < crange_pair[1] and crange_pair[2] < loc_2 and loc_2 < crange_pair[3]):
                        print("\t".join(words))
                        break
            except KeyError:
                continue 
main()
