#!/bin/sh

# align, map


h_brain=DHS_790_hg38
m_brain=DHS_97_mm10


# download gkm-SVM posterior kmer weights
wget https://beerlab.org/gkmalign/human/${h_brain}_300_noproms_nc30_cvAVG_top10k_vs_neg1x_r1_weights_posterior.out -O human_emb_brain_enhancer_weights_posterior.out 
wget https://beerlab.org/gkmalign/mouse/${m_brain}_300_noproms_nc30_cvAVG_top10k_vs_neg1x_r1_weights_posterior.out -O mouse_emb_brain_enhancer_weights_posterior.out
printf "human_emb_brain_enhancer_weights_posterior.out\n" > model_file_emb_brain.txt
printf "mouse_emb_brain_enhancer_weights_posterior.out\n" >> model_file_emb_brain.txt



# run gkm align (x10 faster with precomputed matrix files with -G option)
../../bin/gkm_align -t 1  -g  masker_models.txt -d /mnt/data0/joh27/genomes/ human_mouse_WG_syntenic_intergenic_loci.to_align -W 0.5,model_file_emb_brain.txt -p 50 -o output_files/ -n gkmsvm_emb_brain_weighted_w_0.5

