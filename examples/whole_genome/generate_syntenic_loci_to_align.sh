#!/bin/sh

# 10/25/2023
# Outputs of this pipeline will soon be uploaded to the beerlab/gkm-align webpage


printf "\nThis pipeline first uses LASTZ software to identify short sequence matches in human and mouse. \n"
printf "This step (>5days) can be skipped by downloading a pre-computed intermediate file. \n"
printf "\nPress 1 to skip and download intermediate file from beerlab.org/gkmalign/ (recommended).\n"
printf "Press 2 to run the whole pipeline. LASTZ will be downloaded/compiled and hg38/mm10 multifasta files will be downloaded.\n"
read -p "1 recommended.   " choice
case $choice in
    [1]* )
	    wget https://beerlab.org/gkmalign/short_sequence_human-mouse_syntenic_intergenic.txt -O short_sequence_human-mouse_syntenic_intergenic.txt
            break;;
    [2]* ) 
	    genome_1="hg38"
            genome_2="mm10"

            printf "Setting up LASTZ"
	    wget https://github.com/lastz/lastz/archive/refs/tags/1.04.22.zip -O 1.04.22.zip
	    unzip 1.04.22.zip
	    cd lastz-1.04.22/
	    make
            make lastz_32
	    cd .. 

	    printf "\n\n  downloading hg38 and mm10"
	    wget https://hgdownload.cse.ucsc.edu/goldenpath/${genome_1}/bigZips/${genome_1}.fa.gz -O ${genome_1}.fa.gz
	    wget https://hgdownload.cse.ucsc.edu/goldenpath/${genome_2}/bigZips/${genome_2}.fa.gz -O ${genome_2}.fa.gz
	    gunzip *fa.gz

	    pwd
	    printf "\n\n generating hg38/mm10 short sequence matches. (make take days)\n" 
	    lastz-1.04.22/src/lastz_32 /mnt/data0/joh27/genomes/${genome_2}/multifasta/${genome_2}.fa[multiple] /mnt/data0/joh27/genomes/${genome_1}/multifasta/${genome_1}.fa[multiple] --format=axt --gfextend --nochain --nogapped  --notransition --seed=match10 --step=10 > lastz_output.axt


	    cat lastz_output.axt |grep  chr|grep -v _ | awk -v g1="${genome_1}" -v g2="${genome_2}" '{print g2"\t"$2"\t"$3"\t"$4"\t"g1"\t"$5"\t"$6"\t"$7"\t"$8}' > short_sequence_matches_1.axt

	    python ../../scripts/flip_query_coordinates.py short_sequence_matches_1.axt hg38 /mnt/data0/joh27/projects/alignment_enhancer_conservation/chrom_sizes/hg38.chrom.sizes > short_sequence_matches_2.axt

            # filter using ENCODE syntenic integenic loci
            wget https://beerlab.org/gkmalign/Supplementary_Table_6.txt -O hg38_mm10_ortholog_syntenic_intergenic_pairs.txt
            python ../../scripts/filter_short-seq-matches_by_syntenic_integenic_loci.py  short_sequence_matches_2.axt  hg38_mm10_ortholog_syntenic_intergenic_pairs.txt > short_sequence_human-mouse_syntenic_intergenic.txt

            break ;;
    * ) printf "Please enter  1 or 2\n";;
esac


printf "\n\n\nSetup complete: short sequence matches within hg38/mm10 syntenic intergenic loci\n"
printf "The rest of the pipeline will take less than 20 minutes.\n"

# chain to generate syntenic loci to align
python ../../scripts/chain_short_seq_matches.py  short_sequence_human-mouse_syntenic_intergenic.txt 0.00001 100000 0.1 0 short_sequence_human-mouse_syntenic_intergenic.chains

python ../../scripts/convert_chain_to_to-align.py short_sequence_human-mouse_syntenic_intergenic.chains  20000 10000 1000 human_mouse_WG_syntenic_intergenic_loci.to_align

printf "\n\n Pipeline complete; Saved as: human_mouse_WG_syntenic_intergenic_loci.to_align"
