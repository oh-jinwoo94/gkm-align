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

# Running gkm-align
In this section, we show how we can use gkm-align to run whole-genome alignment between human and mouse. WG alignment between other mammals can also be computed similarly. Running gkm-align requires a pre-computed files containing 1) a list of syntenic intergenic loci and 2) gkm-SVM models for human/mouse genomic background. These files can be found in this repository, and we explain how they can be computed in later sections of this document. 

# Authors
- Jin Woo Oh 
- Michael A. Beer
