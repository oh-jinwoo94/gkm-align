#!/bin/sh


echo "../../data/mouse_genomic_background_model_p_0.1.out" > masker_models.txt
echo "../../data/human_genomic_background_model_p_0.1.out" >> masker_models.txt


mkdir output_files
../../bin/gkm_align  -t 1  HBB.to_align -d ../../data/genomes/ -g masker_models.txt   -p 50 -o output_files -n HBB_LCR_mm10-hg38
../../bin/gkm_align  -t 2  HBB_LCR_enhancers_mm10.bed  -c output_files/HBB_LCR_mm10-hg38.coord -q mm10 -m -o output_files -n HBB_LCR_enhancers_mm10_mapped_to_hg38 



