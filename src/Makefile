# Makefile
CXX = g++
COMMON_FLAGS = -std=c++11 -Wall -Wextra -pedantic -pthread -O3
BIN_DIR = ../bin

# --- OS & Architecture Detection ---
UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

# 1. Detect Hardware Features
# Linux uses lscpu; Mac (Darwin) implies hardware based on arch
HAS_AVX2 := $(shell lscpu 2>/dev/null | grep -q avx2 && echo 1)
HAS_SSE2 := $(shell lscpu 2>/dev/null | grep -q sse2 && echo 1)
IS_ARM64 := $(shell echo $(UNAME_M) | grep -q 'arm64' && echo 1)

# 2. Set Flags
ifeq ($(HAS_AVX2),1)
    # Standard Intel/AMD Server
    CXXFLAGS = $(COMMON_FLAGS) -mavx2
    SIMD_DEFINE = -D__AVX2__
    MSG = AVX2 (Intel/AMD)
else ifeq ($(HAS_SSE2),1)
    # Older Intel
    CXXFLAGS = $(COMMON_FLAGS) -msse2
    SIMD_DEFINE = -D__SSE2__
    MSG = SSE2 (Intel Fallback)
else ifeq ($(IS_ARM64),1)
    # Apple Silicon (M1/M2/M3) or AWS Graviton
    CXXFLAGS = $(COMMON_FLAGS)
    SIMD_DEFINE = -D__ARM_NEON
    MSG = ARM NEON (Apple Silicon)
else
    # Ancient Hardware
    CXXFLAGS = $(COMMON_FLAGS)
    SIMD_DEFINE = -D__NOSIMD__
    MSG = No SIMD
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
