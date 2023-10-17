(will be updated soon, before 10/20/2023, with more details and reproducible examples)
# gkm-align
gkm-align is a whole-genome alignment algorithm specifically designed to map distal enhancers conserved between distant mammals (e.g., human and mouse). gkm-align discovers orthologous enhancers by identifying optimal alignment paths with maximal similarities in gapped-kmer compositions along syntenic loci. gkm-align's performance can further be boosted by incorporating conserved enhancer vocabularies obtained using gkm-SVM sequence models trained on enhancers. 

Please cite the following paper if you use gkm-align:
Oh JW and Beer MA. Gapped-kmer sequence modeling robustly identifies regulatory vocabularies and distal enhancers conserved between evolutionarily distant mammals. 
[bioRxiv](https://www.biorxiv.org/content/10.1101/2023.10.06.561128v1)

# Installation
gkm-align has been compiled and tested for Linux based operating systems (such as Red Hat, Centos, and Rocky Linux, etc.)

First, download the source code using the following command line:
<pre>
git clone https://github.com/oh-jinwoo94/gkm-align.git
</pre>

Then, compile and set up gkm-align by typing:
<pre>
  chmod +x setup.sh
  ./setup.sh
</pre>

# Running gkm-align whole-genome alignment
In this section, we demonstrate how to use gkm-align to generate whole-genome (WG) alignment between human and mouse (option -t 1). WG alignment between other mammals can also be computed similarly. Running gkm-align requires input files containing **1)** the human and mouse genome  **2)** a list of human/mouse syntenic intergenic loci and **3)** gkm-SVM models for human/mouse genomic background. These files can be found in this repository, and we provide information on how they can be computed in later sections of this document. 

<pre>
cd example
../bin/gkm_align -t 1 -d genomes/ -g ../hg38_mm10/genomic_background_models.txt  ../hg38_mm10/syntenic_loci.to_align -o ofiles/ -n unweighted_hg38_mm10
</pre>
The above command line generates the WG alignment output file: unweighted_hg38_mm10.coord.

- **1)** genomes/ is a directory containing subdirectories hg38/ and mm10/, each containing .fa files for each of the chromosomes. These files can be downloaded using the instructions provided in:

  - https://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/
  - https://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/

- **2)** genomic_background_models.txt contains two lines, each specifying a file name for the human or mouse genome background model.

- **3)** syntenic_loci.to_align contains the list of human/mouse syntenic loci. gkm-align computes a gkm-similarity matrix (G) and the optimal alignment path along the matrix for each line of the file.

Option -G can be added to save gkm-simialritiy matrices in a local directory specified by the -o option. gkm-align will automatially look for relevant pre-computed matrix G's in the local directory unless option -O is provided. This substantially reduces the computation time required for cell-specific gkm-SVM weighted alignment (next section).

Details on other software options can be found by typing:
<pre>
../bin/gkm_align -h
</pre>
For example, adding -p 10 option runs gkm-align with 10 parallel multithreads. 

# Incorporating gkm-SVM enhancer vocabularies into gkm-align 
To run gkm-SVM weighted whole-genome alignment, add -W option followed by a chosen magnitude of cell-specific weighting ("c" in the manuscript Fig.4A, ranging from 0 to 1) and the name of a file containing file names for human and mouse gkm-SVM enhancer models. 

<pre>
../bin/gkm_align -t 1  -g  ../hg38_mm10/genomic_background_models.txt -d genomes/ ../hg38_mm10/syntenic_loci.to_align -W 0.5,gkmSVM_human_mouse_brain_models.txt -o ofiles/ -n brain_weighted_hg38_mm10
</pre>  
The above command line generates the WG alignment output file: brain_weighted_hg38_mm10.coord.

# Mapping human enhancers to the mouse genome.
The -t 1 option generates a WG alignment output coordinate file (e.g., unweighted_hg38_mm10.coord), and -t 2 is used to map human enhancers (e.g., those listed in human_brain_enhancers.bed) to the mouse genome (and vice versa) using the coordinate file (like LiftOver).
<pre>
../bin/gkm_align -t 2 -c unweighted_hg38_mm10.coord  human_brain_enhancers_sub.bed  -o ofiles/ -q hg38  -m -n human_brain_enhancers_mapped_to_mm10
</pre>
-q hg38 specifies that the query enhancers are from the human genome (hg38). Adding -m allows mapping each enhancers to multiple loci in the target genome (equivalent to -multiple for LiftOver).  




# Authors
- Jin Woo Oh *
- Michael A. Beer *
