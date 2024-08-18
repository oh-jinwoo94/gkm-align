(being updated with more detail and examples)

# Table of Contents
- [Introduction](#introduction)
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Running gkm-align](#running-gkm-align)
  - [Local sequence alignment & mapping](#local-sequence-alignment-and-mapping)
    - [example: HBB Locus Control Region](#example-hbb-locus-control-region)
    - [example: FADS gene cluster loci](#example-fads-gene-cluster-loci)
  - [Whole-genome alignment & mapping](#whole-genome-alignment-and-mapping)
    - [Genome alignment](#genome-alignment)
        - [Pre-processing for genome alignment](#pre-processing-for-genome-alignment)
        - [cell-type-independent unweighted genome alignment](#cell-type-independent-unweighted-genome-alignment)
        - [cell-type-specific model-weighted genome alignment](#cell-type-specific-model-weighted-genome-alignment)
    - [Genome-wide mapping](#genome-wide-mapping)
# Introduction
gkm-align is a whole-genome alignment algorithm designed to map distal enhancers conserved between distant mammals (e.g., human and mouse). gkm-align discovers orthologous enhancers by identifying alignment paths with maximal similarity in gapped-kmer compositions along syntenic loci. gkm-align's performance can further be enhanced by incorporating conserved enhancer vocabularies obtained using gkm-SVM sequence models trained on enhancers. 

Please cite the following paper if you use gkm-align:

**Oh, J.W., Beer, M.A.** Gapped-kmer sequence modeling robustly identifies regulatory vocabularies and distal enhancers conserved between evolutionarily distant mammals. **Nature Communications** 15, 6464 (**2024**). https://doi.org/10.1038/s41467-024-50708-z

Also, visit the [gkm-align webpage](https://beerlab.org/gkmalign/) to find useful resource files for running gkm-align. 

# System Requirements
gkm-align software is built for Linux-based operating systems (such as Red Hat, CentOS, and Rocky Linux, etc.).
The package has been tested on the following system:
* Rocky Linux release 8.8 (Green Obsidian).

gkm-align utilizes SIMD parallel computation and requires AVX2 support (to check availability, use: 'lscpu | grep avx2'). However, SIMD is only employed for sequence alignment. Therefore, the software can still be used without AVX2 if you plan to use gkm-align with our precomputed genome alignment output files (e.g., [hg38-mm10_unweighted.coord](https://beerlab.org/gkmalign/hg38-mm10_unweighted.coord)) for mapping conserved enhancers. For example, the -t 1 option requires AVX2, but the -t 2 option can be used without it. 
 
# Installation
First, download the source code using the following command line:
<pre>
git clone https://github.com/oh-jinwoo94/gkm-align.git
</pre>

Then, compile and set up gkm-align using the following command:
<pre>
bash setup.sh
</pre>
The script **1)** compiles gkm-align and **2)** downloads gkm-SVM genomic background models (human & mouse) to the 'data/' directory.

Additionally, if you press y (recommended for following the tutorial more easily), **3)** hg38 and mm10 genomes will be downloaded to the 'data/' directory (approximately 6 gigabytes).

The entire process takes less than 3 minutes. 

# Running  gkm-align
## Local sequence alignment and mapping 

### example: HBB Locus Control Region
In this section, we use gkm-align to align the human and mouse HBB Locus Control Region (HBB-LCR) and map mouse HBB-LCR enhancers to human genome (Oh and Beer, **Figure 3G**). 

Enter the following commands.
<pre>
cd examples/HBB_LCR
bash run_gkmalign.sh
</pre>


The 'run_gkmalign.sh' script contains three parts.

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
  * **'-t 1 HBB.to_align'**: Specifies that gkm-align is in "align" mode and uses the input file HBB.to_align, which contains the genomic coordinates of the HBB Locus Control Regions (LCRs) for both human and mouse.
  * **'-d ../../data/genomes'**: Specifiies directory containing the genome data files for human (hg38) and mouse (mm10). These directories should contain chromosome sequence files (e.g., chr1.fa, chr2.fa).
  * **'-g masker_models.txt'**: Specifies the genome background model to use for repeat masking, which helps reduce alignment errors by masking repetitive sequences.
  * **'-p 50'**: Uses 50 parallel threads to speed up processing. This can be adjusted based on available computational resources.
  * **'-o' and '-n'**: Specify the output directory and output file prefix, respectively.

This step generates 'output_files/HBB_LCR_mm10-hg38.coord', which is used as an input for the following step. 
  
**3)** **Mapping** mouse HBB-LCR enhancers to human. 
<pre>
../../bin/gkm_align  -t 2 -c output_files/HBB_LCR_mm10-hg38.coord  HBB_LCR_enhancers_mm10.bed -q mm10 -m -o output_files -n HBB_LCR_enhancers_mm10_mapped_to_hg3
</pre>
  * **'-t 2 -c output_files/HBB_LCR_mm10-hg38.coord'**: Specifies that gkm-align is in "mapping" mode and uses the output from the alignment step (-t 1) as the coordinate mapping file.
  * **'HBB_LCR_enhancers_mm10.bed'**: Input file containing the mm10 coordinates of mouse HBB-LCR enhancers.
  * **'-c output_files/HBB_LCR_mm10-hg38.coord'**: Uses the output from the alignment step (-t 1) as the coordinate mapping file.
  * **'-q mm10'** : Indicates that the query enhancers.bed file is in the mm10 (mouse) genome.
  * **'-m'**: Allows multiple mappings for enhancers that may align to several locations.
  * **'-o' and '-n'**: Specify the output directory and output file prefix, respectively.



Details on other software options can be found by typing:
<pre>
../../bin/gkm_align -h
</pre>

This process takes between 10 seconds and a few minutes depending on your hardware availability.

### example: FADS gene cluster loci
In this section, we align the human and mouse FADS gene cluster loci. 

Enter the following commands.
<pre>
cd examples/FADS_cluster
bash run_gkmalign.sh
</pre>


'run_gkmalign.sh' in this example is almost identical to the version in the previous HBB-LCR example. 
For example, the 'run_gkmalign.sh' script contains the following command lnes:
<pre>
../../bin/gkm_align  -t 1  FADS_loci.to_align -d ../../data/genomes/ -g masker_models.txt   -p 50 -o output_files -n FADS_loci_mm10-hg38 -G
</pre>
 * Adding '-G' option outputs matrix G (binary) to an output directory specified by -o, for each line in the input file with '.to_align' suffix. Output matrix G file is named automatically based on the genomic ranges from which the matrix was computed. 

The output binary file can be converted into a tab-separated file by running: 
<pre> 
../../bin/binary_matrix_2_tsv output_files/mm10-chr19-10194014-10214169-hg38-chr11-61782802-61802911-diff_strand.matrixG output_files/mm10-chr19-10194014-10214169-hg38-chr11-61782802-61802911-diff_strand.tsv
</pre>

Two visualize FADS locus gkm-align output, run 
<pre>
Rscript ../../scripts/visualize_matrix.R output_files/mm10-chr19-10194014-10214169-hg38-chr11-61782802-61802911-diff_strand.tsv invert
Rscript ../../scripts/visualize_coordinates.R output_files/FADS_loci_mm10-hg38.coord
</pre>

| FADS locus gkmsim matrix (G) | FADS gkm-alignment coordinates |
| ------- | ------- |
| ![FADS locus gkmsim matrix (G)](examples/FADS_cluster/png/fads_matrix.png) | ![FADS gkm-alignment](examples/FADS_cluster/png/fads_coords.png) |

This process takes between 10 seconds and a few minutes depending on your hardware availability.

## Whole-genome alignment and mapping

### Genome alignment
#### Pre-processing for genome alignment
The previous two examples (HBB LCR and FADS loci) demonstrated how gkm-align can align a pair of human and mouse loci when their genomic coordinate ranges are well defined, as below:
<pre>
[HBB LCR]
mm10 chr7 103851395 103883181 hg38 chr11 5267522 5302220 same_strand

[FADS locus]
mm10 chr19 10194014 10214169 hg38 chr11 61782802 61802911 diff_strand*
(*note: 'diff_strand' indicates that the human and mouse loci are inverted relative to each another in their respective genome builds.)
</pre>

gkm-align can be applied at the whole genome level by providing a pre-computed list of conserved syntenic loci. It may consist of a list of flanking windows around known conserved transcription start sites, or it may include a list of predicted syntenic intergenic loci derived from a comprehensive list of short sequence matches, as was done for the gkm-align manuscript. In this part of the tutorial, I will describe how to perform whole-genome alignment with gkm-align using syntenic intergenic loci generated using short sequence matches. 

To generate the list of human-mouse syntenic intergenic loci, run following command lines. 
<pre>
cd examples/whole_genome/ 
bash generate_syntenic_loci_to_align.sh
</pre>
The pipeline encoded in the shell script consists of three parts:
1) The pipeline runs the [**LASTZ**](https://github.com/lastz/lastz/) software to generate a comprehensive list of short sequence matches between human and mouse genomes. This step is computationally intensive, and it may take more than 5 days to run depending on hardware availability. For aligning human and mouse (hg38,mm10), this step can be skipped by downloading the output file we have uploaded ([beerlab](https://beerlab.org/gkmalign/short_sequence_human-mouse_syntenic_intergenic.txt)). The pipeline allows you to choose between the two options, and outputs 'short_sequence_human-mouse_syntenic_intergenic.txt'.

2) To generate syntenic blocks, the pipeline runs the **'chain_short_seq_matches.py'** script, which partitions the generated list of short sequence matches to identify and chain together nearby collinear sequence matches in the 2D coordinate space of human and mouse genomes. The script is based on the algorithm described in [Zhang et al. (1994)](https://www.liebertpub.com/doi/10.1089/cmb.1994.1.217), which we have adapted for more intuitive parameterization and simpler usage. This step generates 'short_sequence_human-mouse_syntenic_intergenic.chains'.

3) The last step of the pipeline is to run **'convert_chain_to_to-align.py'**. Syntenic loci, derived from chaining sequence matches in the previous two steps, tend to be very large, making it computationally intensive to compute the gapped-kmer similarity matrices. This code helps expedite the process by breaking the syntenic blocks into smaller pieces through K-means clustering of the sequence matches based on their 2D human-mouse coordinates within each chain. The number of clusters (k) for each chain is automatically determined using the average target block size. The resulting centroids are then used to define the edges of the smaller syntenic blocks. This step generates 'human_mouse_WG_syntenic_intergenic_loci.to_align', which is then used as input for the gkm-align whole genome alignment.  

The following figure shows an example output from running the pipeline described above (GNA12 inversion locus). **Step 1** generates the dots. **Step 2** chains the dots, with each dot colored according to its assigned chain. **Step 3** generates the rectangles, whcih define the boundaries of syntenic blocks which can then be used as input for gkm-align whole-genome alignment. 

![GNA12 syntenic blocks](examples/whole_genome/png/gna12_vis.png) 

#### Cell-type-independent unweighted genome alignment

After generating 'human_mouse_WG_syntenic_intergenic_loci.to_align' through the pipeline described above, whole-genome alignment can now be performed by running:
<pre>
bash run_gkmalign.sh
</pre>
The shell script first generates 'masker_models.txt'. It contains file paths to gkm-SVM genomic background models for human and mouse. These files are downloaded to 'data/' upon setting up gkm-align (bash setup.sh). Then the shell script runs gkm-align on the syntenic blocks using the following command:
<pre>
../../bin/gkm_align  -t 1  human_mouse_WG_syntenic_intergenic_loci.to_align -d ../../data/genomes/ -g masker_models.txt   -p 50 -o output_files -n human_mouse_WG_syntenic_intergenic_loci -G
</pre>

Similar to the FADS and HBB examples, this command line will output a '.coord' file for mapping conserved elements between genome builds. The '-G' option locally saves all the gapped-kmer matrices computed in the directory specified with '-o'. Since computing the gapped-kmer matrices is the most time-consuming part of the algorithm, locally saving the matrices allows the whole-genome alignment to be resumed if the process is interrupted. This option is recommended only if your system has sufficient disk space.

The -G option is particularly useful if you plan to run gkm-align multiple times with various gkm-SVM enhancer models, as it allows the software to reuse previously generated gkm-matrices found in the output directory (-o). For whole genome alignment, this option should be used with discretion because it requires a large disk space (~3TB for hg38-mm10).  Without the locally saved matrices, whole-genome alignment of hg38 and mm10 take about a day with 50 threads ('-p 50'). With the matrices,  whole-genome alignment takes about three hours. 

In general, whole-genome alignment is computationally intensive, and the necessary computational resources may not be available for some users. In such cases, we recommend downloading our pre-computed alignment outputs (analogous to [LiftOver/LASTZ](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/)'s .chain files). These files can be downloaded from our lab website ([beerlab.org](https://beerlab.org/gkmalign/)) under "Output files of aligning the hg38 and mm10 genomes with gkm-align". These files (with suffix .coord) can then be used to map any human elements to mouse or vice versa within conserved intergenic loci, using '-t 2' as described using the [HBB LCR](#example-hbb-locus-control-region). 


#### Cell-type-specific model-weighted genome alignment

### Genome-wide mapping

# Authors
- Jin Woo Oh *
- Michael A. Beer *
