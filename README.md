# gkm-align
gkm-align is a whole-genome alignment algorithm specifically designed to identify distal enhancers conserved between distant mammals such as human and mouse. gkm-align discovers orthologous enhancers by identifying optimal alignment paths with maximal similarities in gapped-kmer compositions along syntenic loci. gkm-align's performance can further be boosted by incorporating conserved enhancer vocabularies obtained using gkm-SVM sequence models trained on enhancers. 

Please cite the following paper if you use gkm-align:
JW Oh, MA Beer  

# Installation
First, download the source code using the folowing command line:
<pre>
git clone https://github.com/oh-jinwoo94/gkm-align.git
</pre>

Then, compile gkm-align by typing:
<pre>
  cd src
  make
</pre>

# Running gkm-align whole-genome alignment
In this section, we demonstrate how to use gkm-align to generate whole-genome (WG) alignment between human and mouse (option -t 1). WG alignment between other mammals can also be computed similarly. Running gkm-align requires input files containing **1)** the human and mouse genome  **2)** a list of human/mouse syntenic intergenic loci and **3)** gkm-SVM models for human/mouse genomic background. These files can be found in this repository, and we provide information on how they can be computed in later sections of this document. 

<pre>
bin/gkm_align -t 1 -d genomes/ -g genomic_background_models.txt  syntenic_loci.2align -o ofiles/ -n unweighted
</pre>
The above command line generates the WG alignment output file: unweighted.coord.

- **1)** genomes/ is a directory containing subdirectories hg38/ and mm10/, each containing .fa files for each of the chromosomes. These files can be downloaded using the instructions provided in:

  - https://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/
  - https://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/

- **2)** genomic_background_models.txt contains two lines, each specifying a file name for the human or mouse genome background model.

- **3)** syntenic_loci.2align contains the list of human/mouse syntenic loci. gkm-align computes a gkm-similarity matrix (G) and the optimal alignment path along the matrix for each line of the file.

Option -G can be added to save gkm-simialritiy matrices in a local directory specified by the -o option. gkm-align will automatially look for relevant pre-computed matrix G's in the local directory unless option -O is provided. This substantially reduces the computation time required for cell-specific gkm-SVM weighted alignment (next section).

Details on other software options can be found by typing:
<pre>
bin/gkm_align -h
</pre>
For example, adding -p 10 option runs gkm-align with 10 parallel multithreads. 

# Running gkm-align whole-genome alignment weighted by gkm-SVM enhancer models. 
To run gkm-SVM weighted whole-genome alignment, add -W option followed by a chosen magnitude of cell-specific weighting ("c" in the manuscript Fig.4A, ranging from 0 to 1) and the name of a file containing file names for human and mouse gkm-SVM enhancer models. 

<pre>
bin/gkm_align -t 1  -g  genomic_background_models.txt -d /mnt/data0/joh27/genomes/ syntenic_loci.2align -W 0.5,gkmSVM_human_mouse_brain_models.txt -o ofiles/ -n brain_weighted
</pre>  
The above command line generates the WG alignment output file: brain_weighted.coord.

# Mapping human enhancers to the mouse genome using the whole-genome alignment output. 
The -t 1 option generates an WG alignment output coordinate file (e.g., unweighted.coord), and -t 2 can be used to map human enhancers to the mouse genome (and vice versa) using the coordinate file (like LiftOver).
<pre>
bin/gkm_align -t 2 -c unweighted.coord human_brain_enhancers.bed  -o ofiles/ -q hg38  -m -n human_brain_enhancers_mapped_to_mm10
</pre>

# Authors
- Jin Woo Oh 
- Michael A. Beer
