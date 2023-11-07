#!/bin/sh


echo "../../data/mouse_genomic_background_model_p_0.1.out" > masker_models.txt
echo "../../data/human_genomic_background_model_p_0.1.out" >> masker_models.txt


mkdir output_files


# align FADS locus. -G option outputs matrix G containing pairwise sequence similarity (gkm-sim matrix)
../../bin/gkm_align  -t 1  FADS_loci.to_align -d ../../data/genomes/ -g masker_models.txt   -p 50 -o output_files -n FADS_loci_mm10-hg38 -G



# For matrix visualization, convert the generated matrix (in binary) to a tab separated file.
../../bin/binary_matrix_2_tsv output_files/mm10-chr19-10194014-10214169-hg38-chr11-61782802-61802911-diff_strand.matrixG output_files/mm10-chr19-10194014-10214169-hg38-chr11-61782802-61802911-diff_strand.tsv

# turn .tsv into .pdf visualization. Output file as ".pdf" appended to the input file name. 
# mm10 FADS locus is inverted relative to hg38 FADS locus. Append "invert" to flip matrix. 
Rscript ../../scripts/visualize_matrix.R output_files/mm10-chr19-10194014-10214169-hg38-chr11-61782802-61802911-diff_strand.tsv invert




# visualize gkm-align coordinates 
Rscript ../../scripts/visualize_coordinates.R output_files/FADS_loci_mm10-hg38.coord
