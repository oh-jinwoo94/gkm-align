#!/bin/sh

to_align=$1 # e.g., human_mouse_WG_syntenic_intergenic_loci.to_align

echo "../../data/mouse_genomic_background_model_p_0.1.out" > masker_models.txt
echo "../../data/human_genomic_background_model_p_0.1.out" >> masker_models.txt

../../bin/gkm_align  -t 1  ${to_align} -d ../../data/genomes/ -g masker_models.txt   -p 50 -o output_files -n hg38-mm10_unweighted -G



