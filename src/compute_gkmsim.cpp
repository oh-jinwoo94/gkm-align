#include "header.h"


// similar functions defined in header, but this version below allows more flexibilities. Defined at the bottom
float compute_matrix(string s1, string s2, int l, vector<unsigned long long int>&  shared_gkm);
float compute_matrix_norm(string s1, string s2, int l, vector<unsigned long long int>& shared_gkm);
unsigned long long int nchoosek(int n, int k);

int main(int argc, char *argv[]) {

	
    if(argc != 5){
        cout << "Usage <0> <ifile name>  <l> <k> <ofname>" <<endl;
        exit(1);
    }

    string ifile_name = argv[1];
    int l = atoi(argv[2]);
    int k = atoi(argv[3]);
    string ofname = argv[4];


    // used for fast gkm sim computation
    vector<unsigned long long int> shared_gkm;
    for (int m = 0; m < l + 1; m++) {
        if ((l - m) >= k) {
            shared_gkm.push_back(nchoosek(l - m, k));
        } else {
            shared_gkm.push_back(0);
        }
    }


    ifstream ifile(ifile_name);

    ofstream ofile;
    if(ofname == "-1"){
	    cout << "output file name must be specified using -o flag" << endl;
	    exit(1);
    }
    ofile.open(ofname);

    cout << "using .. " << endl;
    cout << "l = " << l <<endl;
    cout << "k = " << k << endl;
    string seq1; string seq2;
    string trash;
    string line; string oline;
    while (!ifile.eof()){
        getline(ifile, line);
        if(line.length() == 0){ continue;}
        istringstream ss(line);
        ss >> seq1 >> seq2;
        transform(seq1.begin(), seq1.end(), seq1.begin(), ::toupper);
	transform(seq2.begin(), seq2.end(), seq2.begin(), ::toupper);

        float gkmsim = compute_matrix_norm(seq1, seq2, l, shared_gkm);
	ofile << gkmsim << endl;
    }
    ofile.close();
    ifile.close();
    return 0;
}



float compute_matrix(string s1, string s2, int l, vector<unsigned long long int>&  shared_gkm){
    float tot = 0;
    unsigned int mismatch;

    string::iterator it1 = s1.begin();
    string::iterator it2 = s2.begin();

    for(int i = 0; i <= (int)s1.length()-l; i++){
        for(int j = 0; j <= (int)s2.length()-l; j++){
           mismatch = 0;
           for(int q = 0; q<l; q++){
               char a = *(it1+(i+q));
               char b = *(it2+(j+q));
               mismatch += (a!=b);
           }
           tot += shared_gkm[mismatch];

        }
     }
    return tot;
}




float compute_matrix_norm(string s1, string s2, int l, vector<unsigned long long int>& shared_gkm){
    if((int)s1.length() < l || (int)s2.length() < l){
        return 0;
    }else{
        float k11 = compute_matrix(s1, s1, l, shared_gkm);
        float k22 = compute_matrix(s2, s2, l, shared_gkm);
        float k12 = compute_matrix(s1, s2, l, shared_gkm);
        return k12/(sqrt(k11)*sqrt(k22));
    }
}

unsigned long long int nchoosek(int n, int k) {
    unsigned long long int num = 1;
    unsigned long long int den = 1;

    for (int i = k + 1; i <= n; i++) {
        num = num * i;
    }

    for (int i = 1; i <= (n - k); i++) {
        den = den * i;
    }
    return num / den;
}

