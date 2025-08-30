CC = g++
# Check OS first - only allow Linux
OS_CHECK = $(shell uname -s)
ifeq ($(OS_CHECK),Linux)
    # Check for AVX2 support first, fall back to SSE2, then no SIMD
    AVX2_CFLAGS = -std=c++11 -Wall -Wextra -pedantic -pthread -mavx2
    SSE2_CFLAGS = -std=c++11 -Wall -Wextra -pedantic -pthread -msse2
    NOSIMD_CFLAGS = -std=c++11 -Wall -Wextra -pedantic -pthread

    # Try to compile with AVX2 first, then SSE2, then no SIMD
    AVX2_TEST = $(shell echo 'int main() { __m256i x; return 0; }' | g++ -mavx2 -o /dev/null - 2>/dev/null && echo "AVX2" || echo "SSE2")
    SSE2_TEST = $(shell echo 'int main() { __m128i x; return 0; }' | g++ -msse2 -o /dev/null - 2>/dev/null && echo "SSE2" || echo "NOSIMD")

    ifeq ($(AVX2_TEST),AVX2)
        CFLAGS = $(AVX2_CFLAGS)
        SIMD_DEFINE = -D__AVX2__
    else ifeq ($(SSE2_TEST),SSE2)
        CFLAGS = $(SSE2_CFLAGS)
        SIMD_DEFINE = -D__SSE2__
    else
        CFLAGS = $(NOSIMD_CFLAGS)
        SIMD_DEFINE = -D__NOSIMD__
    endif
else
    $(error This software is designed for Linux systems only. Detected OS: $(OS_CHECK). Please use a Linux system with x86/x64 architecture.)
endif

all: gkm_align  mask_fa compute_gkmsim binary_matrix_2_tsv visualize_mappings train_masker

# Show SIMD support information
	@echo "Compiling with $(AVX2_TEST) support"
	@echo "CFLAGS: $(CFLAGS)"
	@echo "SIMD_DEFINE: $(SIMD_DEFINE)"


gkm_align:  gkm_align.o Matrix.o MatrixG_Computer.o Seq_Aligner.o   Mapper.o functions.o
	$(CC) $(CFLAGS) gkm_align.o Matrix.o  MatrixG_Computer.o Seq_Aligner.o  Mapper.o functions.o -o gkm_align
	mv gkm_align ../bin
gkm_align.o: gkm_align.cpp header.h
	$(CC) $(CFLAGS) $(SIMD_DEFINE) -c gkm_align.cpp

Matrix.o: Matrix.cpp header.h
	$(CC) $(CFLAGS) $(SIMD_DEFINE) -c Matrix.cpp

MatrixG_Computer.o: MatrixG_Computer.cpp header.h
	$(CC) $(CFLAGS) $(SIMD_DEFINE) -c MatrixG_Computer.cpp

Seq_Aligner.o: Seq_Aligner.cpp header.h
	$(CC) $(CFLAGS) $(SIMD_DEFINE) -c Seq_Aligner.cpp

Mapper.o: Mapper.cpp header.h
	$(CC) $(CFLAGS) $(SIMD_DEFINE) -c Mapper.cpp

functions.o: functions.cpp header.h
	$(CC) $(CFLAGS) $(SIMD_DEFINE) -c functions.cpp


mask_fa:  mask_fa.o functions.o
	$(CC) $(CFLAGS) mask_fa.o functions.o -o mask_fa
	mv mask_fa ../bin
mask_fa.o: mask_fa.cpp header.h
	$(CC) $(CFLAGS) $(SIMD_DEFINE) -c mask_fa.cpp

binary_matrix_2_tsv:  binary_matrix_2_tsv.o Matrix.o functions.o
	$(CC) $(CFLAGS) binary_matrix_2_tsv.o functions.o Matrix.o -o binary_matrix_2_tsv
	mv binary_matrix_2_tsv ../bin

binary_matrix_2_tsv.o: binary_matrix_2_tsv.cpp header.h
	$(CC) $(CFLAGS) $(SIMD_DEFINE) -c binary_matrix_2_tsv.cpp

compute_gkmsim:  compute_gkmsim.o
	$(CC) $(CFLAGS) compute_gkmsim.o -o compute_gkmsim
	mv compute_gkmsim ../bin



compute_gkmsim.o: compute_gkmsim.cpp header.h
	$(CC) $(CFLAGS) $(SIMD_DEFINE) -c compute_gkmsim.cpp

visualize_mappings:  visualize_mappings.o Matrix.o MatrixG_Computer.o Seq_Aligner.o  functions.o
	$(CC) $(CFLAGS) visualize_mappings.o  Matrix.o  MatrixG_Computer.o Seq_Aligner.o  functions.o  -o visualize_mappings
	mv visualize_mappings ../bin

visualize_mappings.o: visualize_mappings.cpp header.h
	$(CC) $(CFLAGS) $(SIMD_DEFINE) -c visualize_mappings.cpp




train_masker:  functions.o train_masker.o
	$(CC) $(CFLAGS) train_masker.o functions.o -o train_masker
	mv train_masker ../bin

train_masker.o: train_masker.cpp header.h
	$(CC) $(CFLAGS) $(SIMD_DEFINE) -c train_masker.cpp


clean:
	rm -f *.o gkm_align 

