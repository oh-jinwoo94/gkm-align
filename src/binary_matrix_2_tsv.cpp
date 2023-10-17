#include "header.h"


int main(int argc, char *argv[]) {

	
    if(argc != 3){
        cout << "Usage <0> <input binary matrix file>  <ofname>" <<endl;
        exit(1);
    }

    string ifile_name = argv[1];
    string ofname = argv[2];

    Matrix M = Matrix(ifile_name);
    M.save_matrix_tsv(ofname);

    return 0;
}


