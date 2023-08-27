# gkm-align
gkm-align is a whole-genome alignment algorithm specifically designed to identify distal enhancers conserved between distant mammals such as human and mouse. gkm-align discover orthologous enhancers by identifying optimal alignment paths with maximal similarities in gapped-kmer compositions along syntenic intergenic loci. gkm-align's performance can further be boosted by incorporating conserved enhancer vocabularies obtained using gkm-SVM sequence models trained on enhancers. 

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
In this section, we show how we can use gkm-align to run whole-genome alignment between human and mouse. WG alignment between other mammals can also be computed similarly. Running gkm-align requires input files files containing **1)** the human and mouse genome  **2)** a list of human/mouse syntenic intergenic loci and **3)** gkm-SVM models for human/mouse genomic background. These files can be found in this repository, and we provide information on how they can be computed in later sections of this document. 

<pre>
bin/gkm_align -t 1 -d genomes/ -g genomic_background_models.txt  syntenic_loci.2align -o ofiles/ -n unweighted
</pre>

- **1)** genomes/ is a directory containing subdirectories hg38/ and mm10/, each containing .fa files for each chromosomes. These files can be downloaded using the instructions provided in:

  - https://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/
  - https://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/

- **2)** genomic_background_models.txt contains two lines, each specifying the file name for human and mouse genome background model.

- **3)** syntenic_loci.2align contains a list of human/mouse syntenic loci. gkm-align computes a gkm-similarity matrix (G) and an optimal alignment path along the matrix for each line of this file.
 
# Authors
- Jin Woo Oh 
- Michael A. Beer
