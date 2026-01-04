#include "header.h"


// constructor
Matrix::Matrix(size_t nr, size_t nc): n_rows(nr),  n_cols(nc),  data(nr * nc)
{}

Matrix::Matrix(size_t nr, size_t nc, float val): n_rows(nr),  n_cols(nc),  data(nr * nc, val)
{}

// binary file containing (int)num_row, (int)num_col, val1, val2, ...
Matrix::Matrix(string fname){
	check_file(fname);
	ifstream ifile;
	ifile.open(fname, ios::in | ios::binary);
	ifile.read(reinterpret_cast<char*>(&n_rows), sizeof(size_t));
	ifile.read(reinterpret_cast<char*>(&n_cols), sizeof(size_t));
	data.resize(n_rows*n_cols);

	ifile.read(reinterpret_cast<char*>(&data[0]), n_rows*n_cols*sizeof(float));
	ifile.close();
}


//  update val
float& Matrix::operator()(size_t i, size_t j){
    return data[i * n_cols + j];
}

// return val
float Matrix::operator()(size_t i, size_t j) const{
    return data[i * n_cols + j];
}

vector<size_t> Matrix::dims(){
    return {n_rows, n_cols};
}

float Matrix::compute_mean(){
    float mean = 0;
    for (size_t i = 0; i < n_rows; i++) {
        for (size_t j = 0; j < n_cols; j++) {
            mean = mean + data[i * n_cols + j];
        }
    }
    return mean/(n_cols*n_rows);
}

float Matrix::compute_var(){
    float var = 0;
    float mean = this->compute_mean();
    for (size_t i = 0; i < n_rows; i++) {
        for (size_t j = 0; j < n_cols; j++) {
            var = var + pow(data[i * n_cols + j] - mean, 2);
        }
    }
    var = var / (n_rows * n_cols - 1);
    return var;
}



// save matrix in binary format 
int Matrix::save_matrix(string fname){
	ofstream ofile;   
	ofile.open(fname, ios::out | ios::binary);
	ofile.write(reinterpret_cast<char*>(&n_rows), sizeof(size_t));
	ofile.write(reinterpret_cast<char*>(&n_cols), sizeof(size_t));
	ofile.write(reinterpret_cast<char*>(&data[0]), data.size()*sizeof(float)); 
	ofile.close();
	return 0;
}


// save in human readable .tsv format 
int Matrix::save_matrix_tsv(string fname){
	ofstream ofile;
        ofile.open(fname);
	for(unsigned int i = 0; i < n_rows; i++){
		for(unsigned int j = 0; j < n_cols; j++){
			ofile << data[i*n_cols + j] << "\t";
		}
		ofile << "\n";
	}
        ofile.close();
        return 0;
}





