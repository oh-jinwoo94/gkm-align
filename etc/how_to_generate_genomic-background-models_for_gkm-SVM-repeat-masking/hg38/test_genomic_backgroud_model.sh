#!/bin/sh


# Testing genomic background models with gkm-align by masking fasta files.


# Mask human sample sequences using the human background model (downloaded via setup.py).
../../../bin/mask_fa files/test_hg38.fa ../../../data/human_genomic_background_model_p_0.1.out test_hg38_gkm-masked_h.fa

# Mask human sample sequences using the mouse background model.
../../../bin/mask_fa files/test_hg38.fa ../../../data/mouse_genomic_background_model_p_0.1.out test_hg38_gkm-masked_m.fa


# Comparing the similarity between human and mouse background models.
# gkm-SVM background model identifies overrepresented simple sequence patterns in genomic sequences,
# and most of these patterns are species-invariant. So, if generating these background models is difficult for you,
# you can also consider using the background models we have uploaded (e.g., human-mouse average model for masking another species).

python compare_fa.py  # output: 95.19% match

