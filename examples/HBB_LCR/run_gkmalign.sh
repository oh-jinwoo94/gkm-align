#!/bin/sh


echo "../../data/mouse_genomic_background_model_p_0.1.out" > masker_models.txt
echo "../../data/human_genomic_background_model_p_0.1.out" >> masker_models.txt


mkdir output_files
../../bin/gkm_align  -t 1 -d ../../data/genomes/ -g masker_models.txt   -p 50 -o output_files -n HBB_LCR_mm10-hg38 HBB.to_align
../../bin/gkm_align  -t 2 -c output_files/HBB_LCR_mm10-hg38.coord -q mm10 -m -o output_files -n HBB_LCR_enhancers_mm10_mapped_to_hg38 HBB_LCR_enhancers_mm10.bed



