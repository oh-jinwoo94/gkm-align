#!/bin/sh

# 1. Detect OS and Architecture
OS_CHECK=$(uname -s)
ARCH_CHECK=$(uname -m)

echo "Detected System: $OS_CHECK ($ARCH_CHECK)"

if [ "$OS_CHECK" != "Linux" ] && [ "$OS_CHECK" != "Darwin" ]; then
    echo "Error: This software is designed for Linux or macOS."
    echo "Detected: $OS_CHECK"
    exit 1
fi

# 2. Check Compiler & SIMD Support
echo "Checking compilation environment..."

# Define a temporary test file
TEST_SRC="simd_check.cpp"

if echo "$ARCH_CHECK" | grep -q "arm64"; then
    # --- ARM Logic (Apple Silicon / Graviton) ---
    echo "ARM64 architecture detected. Verifying NEON support..."
    cat <<EOF > $TEST_SRC
#include <arm_neon.h>
int main() { uint8x16_t x = vdupq_n_u8(0); return 0; }
EOF
    # Try compiling (Mac clang defaults to neon, no flags needed usually)
    if ! g++ $TEST_SRC -o /dev/null 2>/dev/null; then
        echo "Warning: Could not compile ARM NEON code."
        echo "Ensure you have Xcode Command Line Tools installed (xcode-select --install)."
        SKIP_ALIGN="true"
    else
        echo "NEON support confirmed."
    fi

else
    # --- Intel/AMD Logic (x86_64) ---
    echo "x86_64 architecture detected. Verifying SSE2/AVX2 support..."
    cat <<EOF > $TEST_SRC
#include <emmintrin.h>
int main() { __m128i x = _mm_setzero_si128(); return 0; }
EOF
    # Try compiling with SSE2 flag
    if ! g++ -msse2 $TEST_SRC -o /dev/null 2>/dev/null; then
        echo "Warning: This system does not support SSE2 SIMD instructions."
        SKIP_ALIGN="true"
    else
        echo "SSE2/AVX2 support confirmed."
    fi
fi

# Clean up test file
rm -f $TEST_SRC

if [ "$SKIP_ALIGN" = "true" ]; then
    echo "gkm-align will compile, but Alignment Mode (-t 1) will be disabled due to lack of SIMD."
    echo "Mapping Mode (-t 2) will still work."
fi

echo "Proceeding with compilation..."

# 3. Compile
cd src
make clean
make
cd ..

# --- Helper Function for Downloads (Handles wget vs curl) ---
download_file() {
    # Usage: download_file <URL> <OUTPUT_NAME> (Optional)
    URL=$1
    OUT=$2
    
    if command -v wget >/dev/null 2>&1; then
        if [ -z "$OUT" ]; then
            wget -q "$URL"
        else
            wget -q "$URL" -O "$OUT"
        fi
    elif command -v curl >/dev/null 2>&1; then
        if [ -z "$OUT" ]; then
            curl -s -L -O "$URL"
        else
            curl -s -L "$URL" -o "$OUT"
        fi
    else
        echo "Error: Neither wget nor curl found. Cannot download files."
        exit 1
    fi
}

# 4. Download Background Models
mkdir -p data
printf "\nDownloading hg38 and mm10 genome background models to data/...\n"

download_file "https://beerlab.org/gkmalign/masker_model_hg38_outside_union_repr_kmer_cluster_threshold_pieces_p_0.1.out" "data/human_genomic_background_model_p_0.1.out"
download_file "https://beerlab.org/gkmalign/masker_model_mm10_outside_union_repr_kmer_cluster_threshold_pieces_p_0.1.out" "data/mouse_genomic_background_model_p_0.1.out"

echo "Models downloaded."

# 5. Download Genomes (Optional Loop)
while true; do
    printf "\nDownload hg38 and mm10 genomes? (y/n) \n"
    read -p "y recommended for the github tutorial (~6 gigabytes).   " yn
    case $yn in
        [Yy]* )
            mkdir -p data/genomes/hg38
            mkdir -p data/genomes/mm10

            # Download hg38
            cd data/genomes/hg38
            printf "\nDownloading hg38 to data/genomes/hg38/...\n"
            for c in {1..22} X Y M 
            do
                printf "chr${c} "
                download_file "https://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr${c}.fa.gz"
            done
            printf "\nUnzipping hg38...\n"
            gunzip -f *.fa.gz

            # Download mm10
            cd ../mm10
            printf "\nDownloading mm10 to data/genomes/mm10/...\n"
            for c in {1..19} X Y M
            do
                printf "chr${c} "
                download_file "https://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chr${c}.fa.gz"
            done
            printf "\nUnzipping mm10...\n"
            gunzip -f *.fa.gz
            
            cd ../../../
            printf "\nDownload complete.\n"
            break;;
            
        [Nn]* ) 
            exit;;
            
        * ) 
            printf "Please enter y or n.\n";;
    esac
done
