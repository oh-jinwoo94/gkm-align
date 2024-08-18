#!/bin/sh

# dowload human brain enhancers
wget https://beerlab.org/gkmalign/human_enh/DHS_790_hg38_300_noproms_nc30.bed

# download the hg38_mm10 alignmnet output file  (can also be generated using run_gkmalign.sh)
wget https://beerlab.org/gkmalign/hg38-mm10_unweighted.coord


# format query enhancers
awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' DHS_790_hg38_300_noproms_nc30.bed > human_brain_enhancers.bed

# run mapper 
../../bin/gkm_align  -t 2  human_brain_enhancers.bed  -c hg38-mm10_unweighted.coord -q mm10 -m -o output_files -n human_brain_enhancers_mapped_to_mm10 
