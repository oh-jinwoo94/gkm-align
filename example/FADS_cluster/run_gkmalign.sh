#!/bin/sh


echo "../../data/mouse_genomic_background_model_p_0.1.out" > masker_models.txt
echo "../../data/human_genomic_background_model_p_0.1.out" >> masker_models.txt


mkdir output_files

# align FADS locus. -G option outputs matrix G containing pairwise sequence similarity
../../bin/gkm_align  -t 1  FADS_loci.to_align -d ../../data/genomes/ -g masker_models.txt   -p 50 -o output_files -n FADS_loci_mm10-hg38 -G

# convert the generated matrix (in binary) to a tab separated file.
../../bin/binary_matrix_2_tsv output_files/mm10-chr19-10194014-10214169-hg38-chr11-61782802-61802911-diff_strand.matrixG output_files/mm10-chr19-10194014-10214169-hg38-chr11-61782802-61802911-diff_strand.tsv

