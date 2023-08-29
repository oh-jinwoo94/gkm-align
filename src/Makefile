CC = g++
CFLAGS = -std=c++11 -Wall -Wextra -pedantic -pthread -mavx2 

all: gkm_align  mask_fa compute_gkmsim binary_matrix_2_tsv visualize_mappings


gkm_align:  gkm_align.o Matrix.o MatrixG_Computer.o Seq_Aligner.o   Mapper.o functions.o
	$(CC) $(CFLAGS) gkm_align.o Matrix.o  MatrixG_Computer.o Seq_Aligner.o  Mapper.o functions.o -o gkm_align
	mv gkm_align ../bin
gkm_align.o: gkm_align.cpp header.h
	$(CC) $(CFLAGS) -c gkm_align.cpp

Matrix.o: Matrix.cpp header.h
	$(CC) $(CFLAGS) -c Matrix.cpp

MatrixG_Computer.o: MatrixG_Computer.cpp header.h
	$(CC) $(CFLAGS) -c MatrixG_Computer.cpp

Seq_Aligner.o: Seq_Aligner.cpp header.h
	$(CC) $(CFLAGS) -c Seq_Aligner.cpp

Mapper.o: Mapper.cpp header.h
	$(CC) $(CFLAGS) -c Mapper.cpp

functions.o: functions.cpp header.h
	$(CC) $(CFLAGS) -c functions.cpp


mask_fa:  mask_fa.o functions.o
	$(CC) $(CFLAGS) mask_fa.o functions.o -o mask_fa
	mv mask_fa ../bin
mask_fa.o: mask_fa.cpp header.h
	$(CC) $(CFLAGS) -c mask_fa.cpp

binary_matrix_2_tsv:  binary_matrix_2_tsv.o Matrix.o functions.o
	$(CC) $(CFLAGS) binary_matrix_2_tsv.o functions.o Matrix.o -o binary_matrix_2_tsv
	mv binary_matrix_2_tsv ../bin

binary_matrix_2_tsv.o: binary_matrix_2_tsv.cpp header.h
	$(CC) $(CFLAGS) -c binary_matrix_2_tsv.cpp

compute_gkmsim:  compute_gkmsim.o
	$(CC) $(CFLAGS) compute_gkmsim.o -o compute_gkmsim
	mv compute_gkmsim ../bin



compute_gkmsim.o: compute_gkmsim.cpp header.h
	$(CC) $(CFLAGS) -c compute_gkmsim.cpp

visualize_mappings:  visualize_mappings.o Matrix.o MatrixG_Computer.o Seq_Aligner.o  functions.o
	$(CC) $(CFLAGS) visualize_mappings.o  Matrix.o  MatrixG_Computer.o Seq_Aligner.o  functions.o  -o visualize_mappings
	mv visualize_mappings ../bin

visualize_mappings.o: visualize_mappings.cpp header.h
	$(CC) $(CFLAGS) -c visualize_mappings.cpp

clean:
	rm -f *.o gkm_align 

