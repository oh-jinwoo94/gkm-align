#!/bin/sh


# specify alignment loci
to_align=$1 # e.g., human_mouse_WG_syntenic_intergenic_loci.to_align

# model specification 
h_model=$2  # e.g., DHS_790_hg38   # for human embryonic brain enhancer model (https://beerlab.org/gkmalign/Supplementary_Table_2.txt)
m_model=$3  # e.g,. DHS_97_mm10    # for mouse emtryonic brain enhancer model (https://beerlab.org/gkmalign/Supplementary_Table_3.txt) 

# enhancer model weighting 
c=$4        # e.g., 0.5            


# download gkm-SVM posterior kmer weights
wget https://beerlab.org/gkmalign/human/${h_model}_300_noproms_nc30_cvAVG_top10k_vs_neg1x_r1_weights_posterior.out -O ${h_model}_enhancer_weights_posterior.out 
wget https://beerlab.org/gkmalign/mouse/${m_model}_300_noproms_nc30_cvAVG_top10k_vs_neg1x_r1_weights_posterior.out -O ${m_model}_enhancer_weights_posterior.out
printf "${h_model}_enhancer_weights_posterior.out\n" > sequence_model_files.txt
printf "${m_model}_enhancer_weights_posterior.out\n" >> sequence_model_files.txt



# run gkm align (x10 faster with precomputed matrix files with -G option)
../../bin/gkm_align -t 1  -g  masker_models.txt -d /mnt/data0/joh27/genomes/ human_mouse_WG_syntenic_intergenic_loci.to_align -W ${c},sequence_model_files.txt -p 50 -o output_files/ -n hg38-mm10_enhancer-model-weighted_${hsamp}-${msamp}_c-${c}.coord

