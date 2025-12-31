# Use CXX for C++ compilers (Standard convention)
CXX = g++

# 1. CRITICAL: Added -O3 for performance
# 2. Changed CFLAGS to CXXFLAGS (Standard for C++)
COMMON_FLAGS = -std=c++11 -Wall -Wextra -pedantic -pthread -O3
BIN_DIR = ../bin

# OS Check
OS_CHECK = $(shell uname -s)
ifneq ($(OS_CHECK),Linux)
    $(error This software is designed for Linux systems only.)
endif

# --- SIMD Detection ---
# Note: This detects the BUILD machine's CPU. 
# If you run this on a cluster node with AVX2, the binary might crash 
# if moved to an older node without AVX2.
HAS_AVX2 := $(shell lscpu 2>/dev/null | grep -q avx2 && echo 1)
HAS_SSE2 := $(shell lscpu 2>/dev/null | grep -q sse2 && echo 1)

ifeq ($(HAS_AVX2),1)
    CXXFLAGS = $(COMMON_FLAGS) -mavx2
    SIMD_DEFINE = -D__AVX2__
    MSG = AVX2 (Fastest)
else ifeq ($(HAS_SSE2),1)
    CXXFLAGS = $(COMMON_FLAGS) -msse2
    SIMD_DEFINE = -D__SSE2__
    MSG = SSE2 (Fallback)
else
    CXXFLAGS = $(COMMON_FLAGS)
    SIMD_DEFINE = -D__NOSIMD__
    MSG = No SIMD (Slowest)
endif

# --- Targets ---

# Define all final binaries with their full paths
TARGETS = $(BIN_DIR)/gkm_align \
          $(BIN_DIR)/mask_fa \
          $(BIN_DIR)/compute_gkmsim \
          $(BIN_DIR)/binary_matrix_2_tsv \
          $(BIN_DIR)/visualize_mappings \
          $(BIN_DIR)/train_masker

# Object files needed for the main aligner
ALIGN_OBJS = gkm_align.o Matrix.o MatrixG_Computer.o Seq_Aligner.o Mapper.o functions.o

all: info directory $(TARGETS)

info:
	@echo "----------------------------------------"
	@echo "Building for Linux with $(MSG) support"
	@echo "CXXFLAGS: $(CXXFLAGS)"
	@echo "----------------------------------------"

# Ensure ../bin exists safely
directory:
	@mkdir -p $(BIN_DIR)

# --- Pattern Rules (Reduces Spaghetti Code) ---

# This one rule handles ALL .cpp -> .o compilation
# Added header.h as a dependency for everything (safest bet for this project)
%.o: %.cpp header.h
	$(CXX) $(CXXFLAGS) $(SIMD_DEFINE) -c $< -o $@

# --- Linking Rules ---

$(BIN_DIR)/gkm_align: $(ALIGN_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/mask_fa: mask_fa.o functions.o
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/binary_matrix_2_tsv: binary_matrix_2_tsv.o Matrix.o functions.o
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/compute_gkmsim: compute_gkmsim.o
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/visualize_mappings: visualize_mappings.o Matrix.o MatrixG_Computer.o Seq_Aligner.o functions.o
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/train_masker: train_masker.o functions.o
	$(CXX) $(CXXFLAGS) $^ -o $@

# Clean now cleans the objects AND the binaries in ../bin
clean:
	rm -f *.o $(TARGETS)
