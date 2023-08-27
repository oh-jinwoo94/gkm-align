# gkm-align
gkm-align is a whole-genome alignment algorithm that is specifically designed to identify distal enhancers conserved between distant mammals such as human and mouse. gkm-align discover orthologous enhancers by identifying optimal alignment paths with maximal similarities in gapped-kmer compositions along syntenic intergenic loci. Its performance can further be boosted by incorporating conserved enhancer vocabularies obtained using gkm-SVM sequence models trained on enhancers. 

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


# Authors
Jin Woo Oh

Michael A. Beer
