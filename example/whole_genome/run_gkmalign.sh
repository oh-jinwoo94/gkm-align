#!/bin/sh

# align, map

echo "../../data/mouse_genomic_background_model_p_0.1.out" > masker_models.txt
echo "../../data/human_genomic_background_model_p_0.1.out" >> masker_models.txt

../../bin/gkm_align  -t 1  human_mouse_WG_syntenic_intergenic_loci.to_align -d ../../data/genomes/ -g masker_models.txt   -p 50 -o output_files -n human_mouse_WG_syntenic_intergenic_loci
#../../bin/gkm_align  -t 2  HBB_LCR_enhancers_mm10.bed  -c output_files/HBB_LCR_mm10-hg38.coord -q mm10 -m -o output_files -n HBB_LCR_enhancers_mm10_mapped_to_hg38


