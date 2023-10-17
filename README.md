(being updated with more details and examples)

# Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Running gkm-align](#running-gkm-align)
  - [example: HBB Locus Control Region](#example-hbb-locus-control-region)
  - [example: human-mouse whole-genome alignment](#example-human-mouse-whole-genome-alignment)
  
# Introduction
gkm-align is a whole-genome alignment algorithm designed to map distal enhancers conserved between distant mammals (e.g., human and mouse). gkm-align discovers orthologous enhancers by identifying alignment paths with maximal similarities in gapped-kmer compositions along syntenic loci. gkm-align's performance can further be boosted by incorporating conserved enhancer vocabularies obtained using gkm-SVM sequence models trained on enhancers. 

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
The script **1)** compiles gkm-align and **2)** downloads gkm-SVM genomic background models (human & mouse) to data/.

Further, if you press y (recommended for following the tutorial more easily), **3)** hg38 and mm10 genomes download to data/.

# Running gkm-align

## example HBB Locus Control Region
In this section, we use gkm-align to align the human,mouse HBB Locus Control Region (HBB-LCR) and map mouse HBB-LCR enhancers to human genome (Oh and Beer, **Figure 3G**). 

Enter the following commands.
<pre>
  cd examples/HBB_LCR
  chmod +x run_gkmalign.sh
  ./run_gkmalign.sh
</pre>


'run_gkmalign.sh' script contains three parts.

**1)** Setting up gkm-align input file (specifying gkm-SVM genomic masker models to be used) and output directory. 
<pre>
  echo "../../data/mouse_genomic_background_model_p_0.1.out" > masker_models.txt
  echo "../../data/human_genomic_background_model_p_0.1.out" >> masker_models.txt
  mkdir output_files
</pre>

**2)** **Aligning** human and mouse HBB-LCRs. 
<pre>
  ../../bin/gkm_align  -t 1  HBB.to_align -d ../../data/genomes/ -g masker_models.txt   -p 50 -o output_files -n HBB_LCR_mm10-hg38
</pre>
  * -t 1 option specifies that gkm-align is in 'align' mode.
  * 'HBB.to_align' contains genomic coordinate ranges of human and mouse HBB Locus Control Regions that gkm-align will perform alignment. 
  * -d ../../data/genomes specifiies directory containing hg38/ and mm10/, each containing chromosome.fa. \
  * -g masker_models.txt specifies gkm-SVM genome background model to use for repeat masking.
  * -p 50: use 50 parallel threads. Set to -p 1 if resource unavailable.
  * -o and -n each specifies output directory and output file prefix. 
  
**3)** **Mapping** mouse HBB-LCR enhancers to human. 
<pre>
  ../../bin/gkm_align  -t 2  HBB_LCR_enhancers_mm10.bed  -c output_files/HBB_LCR_mm10-hg38.coord -q mm10 -m -o output_files -n HBB_LCR_enhancers_mm10_mapped_to_hg3
</pre>
  * -t 2 option specifies that gkm-align is in 'mapping' mode.  
  * 'HBB_LCR_enhancers_mm10.bed' contains mm10 coordinates of mouse HBB-LCR enhancers.
  * -c output_files/HBB_LCR_mm10-hg38.coord: output of -t 1 (align mode) to be use for mapping.
  * -q mm10 specifies that the query enhancers.bed is in mm10.
  * -m: allows multiple mapping.
  * -o and -n each specifies output directory and output file prefix. 

Details on other software options can be found by typing:
<pre>
../bin/gkm_align -h
</pre>

## example human-mouse whole-genome alignment

(description coming soon. all necessary command lines are provided in example/whole_genome)

# Authors
- Jin Woo Oh *
- Michael A. Beer *
