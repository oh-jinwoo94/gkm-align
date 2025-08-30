#!/bin/bash

echo "HBB_LCR Test"

# Check if we're in the right directory
if [ ! -f "../src/Makefile" ]; then
    echo "Error: Run from test/ directory"
    exit 1
fi

# Create data directories
mkdir -p ../data/genomes/hg38
mkdir -p ../data/genomes/mm10

# Download background models
echo "Downloading background models..."
wget -q https://beerlab.org/gkmalign/masker_model_hg38_outside_union_repr_kmer_cluster_threshold_pieces_p_0.1.out -O ../data/human_genomic_background_model_p_0.1.out
wget -q https://beerlab.org/gkmalign/masker_model_mm10_outside_union_repr_kmer_cluster_threshold_pieces_p_0.1.out -O ../data/mouse_genomic_background_model_p_0.1.out

# Download chromosomes
echo "Downloading chromosomes..."
cd ../data/genomes/hg38
wget -q "https://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr11.fa.gz"
gunzip chr11.fa.gz
cd ../mm10
wget -q "https://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chr7.fa.gz"
gunzip chr7.fa.gz
cd ../../..

# Check if gkm_align exists
if [ ! -f "../bin/gkm_align" ]; then
    echo "Error: gkm_align not found. Compile first: cd ../src && make"
    exit 1
fi

# Run test
echo "Running test..."
cd ../examples/HBB_LCR
mkdir -p output_files

# Create masker models file
echo "../../data/mouse_genomic_background_model_p_0.1.out" > masker_models.txt
echo "../../data/human_genomic_background_model_p_0.1.out" >> masker_models.txt

# Run gkm-align
../../bin/gkm_align -t 1 HBB.to_align -d ../../data/genomes/ -g masker_models.txt -p 50 -o output_files -n HBB_LCR_mm10-hg38
../../bin/gkm_align -t 2 HBB_LCR_enhancers_mm10.bed -c output_files/HBB_LCR_mm10-hg38.coord -q mm10 -m -o output_files -n HBB_LCR_enhancers_mm10_mapped_to_hg38

# Compare results
echo "Comparing results..."
if diff expected_output/HBB_LCR_enhancers_mm10_mapped_to_hg38.multiple_mapped output_files/HBB_LCR_enhancers_mm10_mapped_to_hg38.multiple_mapped > /dev/null; then
    echo "SUCCESS: Output matches expected results"
    exit 0
else
    echo "FAILURE: Output does not match expected results"
    diff expected_output/HBB_LCR_enhancers_mm10_mapped_to_hg38.multiple_mapped output_files/HBB_LCR_enhancers_mm10_mapped_to_hg38.multiple_mapped
    exit 1
fi
