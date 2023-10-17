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

void Matrix::weigh_rows_by_vect(vector<float> v){
    if(v.size() != n_rows){
	    cout << "dimension mismatch in weigh_rows_by_vect" << endl;
	    exit(1);
    }
    for (size_t i = 0; i < n_rows; i++) {
        for (size_t j = 0; j < n_cols; j++) {
            data[i * n_cols + j] *= v[i];
        }
    }
}

void Matrix::weigh_cols_by_vect(vector<float> v){
    if(v.size() != n_cols){
            cout << "dimension mismatch in weigh_cols_by_vect" << endl;
            exit(1);
    }
    for (size_t j = 0; j < n_cols; j++) {
        for (size_t i = 0; i < n_rows; i++) {
            data[i * n_cols + j] *= v[j];
        }
    }
}




void Matrix::rowwise_avgdotproduct_matrices(vector<vector<float>> M1, vector<vector<float>> M2, string rel_strand_loc){
    if(M1.size() != n_rows || M2.size() != n_cols){
        cout << "dimension mismatch" << endl;
        exit(1);
    }else{

	float dp;
	unsigned int N;
        vector<float> p1; vector<float> p2; 
        for(unsigned int i = 0 ; i < n_rows; i++){
            p1 = M1[i];
            for(unsigned int j = 0; j < n_cols; j++){
                if(rel_strand_loc == "same_strand"){
                    p2 = M2[j];
                }else{
                    p2 = M2[M2.size()-1-j];
                }

		N = p1.size(); 
		dp = 0;
		for(unsigned int k = 0; k < N; k++){
			dp = dp + p1[k]*p2[k];
		}
                data[i*n_cols + j] = dp / N;
            }
        }
    }
}


// row-wise corr (using pearson corr) between two (max-filtered) matrices
void Matrix::rowwise_pcorr_matrices(vector<vector<float>> M1, vector<vector<float>> M2, float max_filter, string rel_strand_loc){
    if(M1.size() != n_rows || M2.size() != n_cols){
        cout << "dimension mismatch" << endl;
        exit(1);
    }else{

	vector<float> p1; vector<float> p2; float max1; float max2;
        for(unsigned int i = 0 ; i < n_rows; i++){
	    p1 = M1[i];
            for(unsigned int j = 0; j < n_cols; j++){
		if(rel_strand_loc == "same_strand"){
		    p2 = M2[j];
		}else{
		    p2 = M2[M2.size()-1-j];
		}

                max1 = *max_element(begin(p1), end(p1));
                max2 = *max_element(begin(p2), end(p2));

                float filter1 = 0 ; float filter2 = 0;
                if(max1>max_filter){
                    filter1 = 1;
                }
                if(max2>max_filter){
                    filter2 = 1;
                }

                data[i*n_cols + j] = (pcorr(p1, p2) + 1 ) /2 * filter1 * filter2;
            }
        }
    }
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





