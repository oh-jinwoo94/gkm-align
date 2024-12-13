#include "header.h"

int main(int argc, char *argv[]) {

    if(argc != 6){
        cout << "Usage <0> <genome_kmer weights> <range file> <genome loc> <fraction mask> <ofile name>" <<endl;
        exit(1);
    }

    string wfile_name = argv[1];
    string rfile_name = argv[2];
    string dir = argv[3];
    float frac  = stof(argv[4]);

    unordered_map<string, float> kmer_weights;
    load_weights(kmer_weights, wfile_name); 


    cout << "train start" << endl;
    float threshold = train_masker(kmer_weights, rfile_name, dir, frac);


    string ofname = argv[5];
    ofstream ofile(ofname);



    ifstream wfile(wfile_name);   
    ofile << "* threshold: " << threshold << endl;
    ofile << "# Threshold obtained using using:" << endl;
    string line;
    ifstream rfile(rfile_name);
    while (getline(rfile, line)) {
        ofile<<"# " << line<<endl;
    }
    
    

    while (getline(wfile, line)) { 
        ofile<<line<<endl;;
    }


    cout << "output file saved as: " << ofname << endl;
    rfile.close();
    wfile.close();
    ofile.close();
    return 0;
}
