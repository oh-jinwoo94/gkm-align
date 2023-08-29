#!/bin/sh

cp /mnt/VolB/data/joh27/alignment_data_storage/new_dataset/noproms_nc30_posterior/cv_posterior/DHS_790_hg38_300_noproms_nc30_cvAVG_top10k_vs_neg1x_r1_weights_posterior.out human_brain_enhancer_model.out

cp /mnt/VolB/data/joh27/alignment_data_storage/new_dataset/noproms_nc30_posterior/cv_posterior/DHS_97_hg38_300_noproms_nc30_cvAVG_top10k_vs_neg1x_r1_weights_posterior.out mouse_brain_enhancer_model.out

sort -grk 4 /data6/mbeer3/data/ENCODE_DHS_bamp/noproms_nc30/DHS_790_hg38_300_noproms_nc30.bed |head -100 > human_brain_enhancers_sub.bed

printf "human_brain_enhancer_model.out\n" > gkmSVM_human_mouse_brain_models.txt
printf "mouse_brain_enhancer_model.out\n" >> gkmSVM_human_mouse_brain_models.txt
