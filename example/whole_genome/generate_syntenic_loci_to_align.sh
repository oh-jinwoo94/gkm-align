#!/bin/sh

# 10/25/2023
# Outputs of this pipeline will soon be uploaded to the beerlab/gkm-align webpage

genome_1="hg38"
genome_2="mm10"

##### run lastz (~10 days) #####  output will be uploaded soon to the beerlab gkmalign website ##
/mnt/data0/joh27/tools/lastz/lastz/src/lastz_32 /mnt/data0/joh27/genomes/${genome_2}/multifasta/${genome_2}.fa[multiple] /mnt/data0/joh27/genomes/${genome_1}/multifasta/${genome_1}.fa[multiple] --format=axt --gfextend --nochain --nogapped  --notransition --seed=match10 --step=1 > lastz_output.axt


axtfile=lastz_output.axt
######## steps below < 20mins  ########
# process lastz output
cat ${axtfile} |grep  chr|grep -v _ | awk -v g1="${genome_1}" -v g2="${genome_2}" '{print g2"\t"$2"\t"$3"\t"$4"\t"g1"\t"$5"\t"$6"\t"$7"\t"$8}' >> short_sequence_matches_1.axt

python ../../scripts/flip_query_coordinates.py short_sequence_matches_1.axt hg38 /mnt/data0/joh27/projects/alignment_enhancer_conservation/chrom_sizes/hg38.chrom.sizes > short_sequence_matches_2.axt

# filter using ENCODE syntenic integenic loci
wget https://beerlab.org/gkmalign/Supplementary_Table_6.txt -O hg38_mm10_ortholog_syntenic_intergenic_pairs.txt

python ../../scripts/filter_short-seq-matches_by_syntenic_integenic_loci.py  short_sequence_matches_2.axt  hg38_mm10_ortholog_syntenic_intergenic_pairs.txt > short_sequence_human-mouse_syntenic_intergenic.axt


# chain to generate syntenic loci to aign
python ../../scripts/chain_short_seq_matches.py  short_sequence_human-mouse_syntenic_intergenic.axt 0.00001 100000 0.1 0 short_sequence_human-mouse_syntenic_intergenic.chains

python ../../scripts/convert_chain_to_to-align.py short_sequence_human-mouse_syntenic_intergenic.chains  20000 10000 1000 human_mouse_WG_syntenic_intergenic_loci.to_align
