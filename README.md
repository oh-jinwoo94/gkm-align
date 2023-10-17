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
The script compiles gkm-align and downloads gkm-SVM genomic background models (human & mouse) to data/.

Further, if you press y (recommended for following the tutorial more easily), hg38 and mm10 genomes download to data/.

# gkm-align tutorial
In this section, we use gkm-align to align the human,mouse HBB Locus Control Region (HBB-LCR) and map mouse HBB-LCR enhancers to human genome (Oh and Beer, **Figure 3G**). 

<pre>
  cd examples/HBB_LCR
  chmod +x run_gkmalign.sh
  ./run_gkmalign.sh
</pre>

<pre>
  echo "../../data/mouse_genomic_background_model_p_0.1.out" > masker_models.txt
echo "../../data/human_genomic_background_model_p_0.1.out" >> masker_models.txt


mkdir output_files
../../bin/gkm_align  -t 1  HBB.to_align -d ../../data/genomes/ -g masker_models.txt   -p 50 -o output_files -n HBB_LCR_mm10-hg38
../../bin/gkm_align  -t 2  HBB_LCR_enhancers_mm10.bed  -c output_files/HBB_LCR_mm10-hg38.coord -q mm10 -m -o output_files -n HBB_LCR_enhancers_mm10_mapped_to_hg38

</pre>




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
