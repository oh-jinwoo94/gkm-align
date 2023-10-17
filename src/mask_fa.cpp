#include "header.h"
// assumes contiguous nucleotides are given in a single line.
// i.e. for genome.fa files, they need to be preprocessed. 

string mask_seq_lc(string seq, const unordered_map<string, float>& kmer_weights, const float threshold) {

    vector<float> bp_scores = gen_sc_v(seq, kmer_weights);
    for (unsigned int i = 0; i < seq.length(); i++) {
        if (bp_scores[i] > threshold) {
            seq[i] = tolower(seq[i]);
        }
    }
    return seq;
}



int main(int argc, char *argv[]) {
    if(argc != 4){
        cout << "Usage <0> <ifile name>  <masker model> <ofile name>" <<endl;
        exit(1);
    }

   string ifile_name = argv[1];
   check_file(ifile_name);
   ifstream ifile(ifile_name);

   unordered_map<string, float> kmer_weights;
   float threshold = load_weights_threshold(kmer_weights, argv[2]);

   ofstream ofile;
   ofile.open(argv[3]);

    string line; 
    while (!ifile.eof()){
        getline(ifile, line);
	if(line[0] == '>'){
	    ofile << line << endl;
	}else{
	   string seq = line;
	   seq = preprocess_seq(seq);
           seq = mask_seq_lc(seq, kmer_weights, threshold); 
	   ofile << seq << endl;
	}
    }
    ofile.close();
    ifile.close();
    return 0;
}



