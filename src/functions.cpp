// some simple and useful functions 


#include "header.h"
const vector<char> ATGC = {'A', 'T', 'G', 'C'};
const unordered_map<char,char> comp = {{'A', 'T'}, {'T', 'A'}, {'G', 'C'}, {'C', 'G'}};

void help(){
	    cout << endl <<endl <<endl;;
            cout << " ========================================================================================================" <<endl;
            cout << "| gkm-Align: gapped-kmer based whole-genome alignment software for mapping conserved distal enhancers     |" <<endl;
            cout << "|                                                                                                         |" << endl;
            cout << "| version: v.1.0                                                                                          |" << endl;
            cout << "|                                                                                                         |" << endl;
            cout << "| For more details, please refer to the github README and the gkm-Align paper.                            |" <<endl;
	    cout << "| - github: https://github.com/oh-jinwoo94/gkm-align                                                      |" << endl;
	    cout << "| - paper :                                                                                               |" << endl;
            cout << " ========================================================================================================\n" << endl;
            cout << "Usage: " <<endl;
            cout << "                  (1) for genome alignment (-t 1):" << endl;
	    cout << "                          ./gkm_align -t 1 -d <genomes> -g <genome_background_models.txt> <input: loci.2align>\n" << endl;
            cout << "                  (2) for enhancer mapping (-t 2):" << endl;
	    cout << "                          ./gkm_align -t 2 -c <alignment_output.coord> -m <input: enhancers.bed> \n" << endl;


            cout << "Software arguments & options" << endl;
            cout << "(1) for genome alignment (-t 1): " << endl;
	    cout << "  input argument      [required] 2align containing interspecies syntenic loci to align." << endl;
            cout << "  -d <string>         [required] directory containing hg38/ and mm10/" <<endl;
            cout << "                                 (or other genomes of interest containing .fasta files for each chromosomes)" << endl;
            cout << "  -g <string>         [required] .txt containing file names for a genome background model for each species\n" << endl;
            cout << "  -l <int>            Gapped-kmer width (default: 11; e.g., -TGA-TCAT--)" << endl;
            cout << "  -k <int>            Number of non-gapped positions in the gapped kmers (default: 7; e.g., -TGA-TCAT--)" << endl;
            cout << "                      - must be smaller than or equal to l." << endl;
            cout << "  -w <int>            Width of the sliding windows (default: 300)" << endl;
            cout << "  -s <int>            Sliding step size of the sliding windows (default: 20)" << endl;
            cout << "                      - w must be divisiable by s (e.g., w = 300 and s = 20)" << endl;
            cout << "  -i <float>          Indel penalty for insertion/deletion (default: 0; scale: Z-score)" << endl;
            cout << "  -p <int>            Number of parallel threads for multithreaded alignment. (default: 1)" <<endl;
	    cout << "  -W <float,string>   Use this option for cell-specific genome-alignment weighted by gkm-SVM enhancer models." << endl;
            cout << "                      float : enhancer model weight ('c': 0-1)" << endl;
            cout << "                      string: .txt containing file names of the enhancer models" <<endl;  
	    cout << "  -G                  Save gkm-similarity matrices to a local directory (default: no save; -o dir/)" << endl;
	    cout << "  -O                  By default, gkm-Align uses relevant gkm-matrices in a directory specified by the -o option." << endl;
	    cout << "                      Provide -O to ignore existing gkm-matrices\n" << endl;

            cout << "(2) for enhancer mapping (-t 2): " << endl;
            cout << "  input argument      [required] .bed containing query enhancers to map." << endl;
            cout << "  -c <string>         [required] .coord output file from -t 1 genome alignment" << endl;
            cout << "  -q <string>         [required] query genome (e.g., hg38 for mapping human enhancers to the mouse genome) \n" << endl;
            cout << "  -m or -u            [required] either -m or -u must be provided. " << endl;
	    cout << "                                 -m: allows mapping enhancers to multiple elements." <<endl;
	    cout << "                                 -u: restrict to unique mapping. \n" << endl; 

	    cout << "shared options: " << endl;
	    cout << "   -o                 output directory   (default: current directory)" <<endl;
	    cout << "   -n                 output file prefix (default: input file prefix )" << endl;

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


void check_file(string fname)
{
    ifstream ifile(fname);
    if(ifile.good() == false){
        cout<< "The following file is not available: " << fname <<endl;
        exit(1);
    }
    ifile.close();
}





// extract directory name from file name. 
// returns "" if no dir name
string dirname(string fname){
    int loc = -1;
    for (unsigned int i = 0; i < fname.length(); i++) {
        if (fname[i] == '/') {
            loc = i;
        }
    }

    string dir = fname.substr(0, loc + 1);
    return dir;
}

string basename(string fname){
    int loc = -1;
    for (unsigned int i = 0; i < fname.length(); i++) {
        if (fname[i] == '/') {
            loc = i;
        }
    }

    string base = fname.substr(loc + 1, fname.length() - (loc+1));
    return base;
}




string getexepath(){
  char result[ PATH_MAX ];
  ssize_t count = readlink( "/proc/self/exe", result, PATH_MAX );
  return dirname(string( result, (count > 0) ? count : 0 ));
}


// if unknown character is given (non ATGC), leave as it is. 
string revcomp(string seq) {
    string oseq = "";
    for (char c : seq) {
        auto it = comp.find(c);
        if (it != comp.end()){
            oseq = (it->second) + oseq;
        }else{
            oseq = c + oseq;
        }
    }
    return oseq;
}

// just turn into upper case and replace unknonw characters with random ATGC
string preprocess_seq(string seq) {
    transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    filter_unknown(seq);
    return seq;
}



// read the entire sequence in the single fasta file
string read_fa(string fname) {
    check_file(fname);
    ifstream fa_file(fname);
    string seq = "";
    string line;
    while (getline(fa_file, line)) {
        seq += line;
    }
    fa_file.close();
    return seq;
}

// read fasta file to extract sequences in the predefined range
// [start, end)
// if start <0, read from beginning of the file. if end > chromosome length, then read to the end 
string read_fa(string fname, int start, int end) {
    check_file(fname);
    ifstream fa_file(fname);
    string seq = "";
    int curr_coord = 0;
    string line;
    bool read_seqs;
    if(start <= 0){
        read_seqs = true;
        start = 0;
    }else{
        read_seqs = false;
    }

    while (getline(fa_file, line)) {
        if (read_seqs == true) {
            if (curr_coord + line.length() > (unsigned int)end) {
                seq += line.substr(0, end - curr_coord);
                break;
            } else {
                seq += line;
            }
        } else if (curr_coord + line.length() > (unsigned int)start) {
            read_seqs = true;
            int s = start - curr_coord;
            seq += line.substr(s, line.length() - s);
        }
        curr_coord += line.length();
    }
    fa_file.close();
    return seq;
}

// assuming all in upper case.
void filter_unknown(string &seq) {
    vector<char>::const_iterator it;
    for (char &c : seq) {
        it = std::find(ATGC.begin(), ATGC.end(), c);
        if (it == ATGC.end()) {
            c = ATGC[rand() % 4];
        }
    }
}



pair<string, string> generate_regions_to_align(string gdir, string build1, string build1_chr, int build1_start, int build1_end,
		                               string build2, string build2_chr, int build2_start, int build2_end,
					       int window, int slide){
    string build1_fname = gdir + build1 + "/" + build1_chr + ".fa";
    string build2_fname = gdir + build2 + "/" + build2_chr + ".fa";

    string seq1 = read_fa(build1_fname, build1_start, build1_end);
    string seq2 = read_fa(build2_fname, build2_start, build2_end);


    seq1 = preprocess_seq(seq1);
    seq2 = preprocess_seq(seq2);


    if ((window % slide) != 0) {
        cout << "window width must be divisible by the slide width" << endl;
        exit(1);
    } // part 1 complete

    if ((seq1.length() - window) % slide != 0) {
        int length_final = window + ((seq1.length() - window) / slide) * slide;
        seq1 = seq1.substr(0, length_final);
    }

    if ((seq2.length() - window) % slide != 0) {
        int length_final = window + ((seq2.length() - window) / slide) * slide;
        seq2 = seq2.substr(0, length_final);
    }

    return make_pair(seq1, seq2);
}

// regular weight file
void load_weights(unordered_map<string, float> &kmer_weights,
                  string fname) {

    check_file(fname);
    kmer_weights.clear();
    ifstream wfile(fname);
    string kmer;
    float weight;
    string line;
    while (getline(wfile, line)) {
        std::istringstream ss(line);
        ss >> kmer >> weight;
        kmer_weights[kmer] = weight;
    }
    wfile.close();
}

tuple<int, float, float, float, float, unsigned int> load_post_weights(unordered_map<string, float> &kmer_weights,
                  string fname) {
        // format:
	// w	300
        // mu_pos  -0.4698601590268965
        // var pos 0.0002796103768422355
       // mu_neg  -0.5784510419876551
        // var_neg 0.0010991364946907067
        // CCCCGATATAC     -0.6253484

    check_file(fname);
    kmer_weights.clear();
    ifstream wfile(fname);
    cout << fname << endl;
    string kmer;
    float weight;
    string line;
    int pred_width; float mu_pos; float var_pos; float mu_neg; float var_neg;
    unsigned int k;
    string trash;
    int i = 0;
    while (getline(wfile, line)) {
	std::istringstream ss(line);
	if(i==0){ss >> trash >> pred_width;}
	if(i == 1){ss >> trash >> mu_pos;}
	else if(i == 2){ss >> trash >> var_pos;}
	else if(i == 3){ss >> trash >> mu_neg;}
	else if(i == 4){ss >> trash >> var_neg;}
	else{
            ss >> kmer >> weight;
            kmer_weights[kmer] = weight;
	    if(i==5){k = kmer.size();}
	}
	i++;
    }
    wfile.close();
    return make_tuple(pred_width, mu_pos, var_pos, mu_neg, var_neg, k);
}


tuple<int, float, float, float, float, unsigned int> load_post_weights_to_matrix(unordered_map<string, vector<float>> &kmer_weight_matrix,
                  string fname) {
//w       300.0
//mu_pos  0.3888576375783236
//var_pos 0.0326371041453618
//mu_neg  0.040440167715889916
//var_neg 0.04457310205165671
//AAAAAAAAAAA     1.5321328531827556


    check_file(fname);
    ifstream wfile(fname);
    cout << fname << endl;
    string kmer;
    float weight;
    string line;
    int pred_width; float mu_pos; float var_pos; float mu_neg; float var_neg; 
    int k;
    string trash;
    int i = 0;
    while (getline(wfile, line)) {
        std::istringstream ss(line);
	if(i==0){ss >> trash >> pred_width;}
        if(i == 1){ss >> trash >> mu_pos;}
        else if(i == 2){ss >> trash >> var_pos;}
        else if(i == 3){ss >> trash >> mu_neg;}
        else if(i == 4){ss >> trash >> var_neg;}
        else{
            ss >> kmer >> weight;
	    if(kmer_weight_matrix.find(kmer) == kmer_weight_matrix.end()){
		kmer_weight_matrix[kmer] = {weight};
	    }else{
                kmer_weight_matrix[kmer].push_back(weight);
	    }
            if(i==5){k = kmer.size();}
        }
        i++;
    }
    wfile.close();
    return make_tuple(pred_width,  mu_pos, var_pos, mu_neg, var_neg, k);
}

// this is for Masker. has both kmer weights and a threshold value. 
// format: 
// * threshold: 0.359557
// # Threshold obtained using using:
// # hg38    chr21   20000000 20100000
// AAAAAAAAAAA     1.6789377
// ATATATATATA     1.5665055
// ACACACACACA     1.208083756
float load_weights_threshold(unordered_map<string, float> &kmer_weights,
                  string fname) {

    check_file(fname);
    kmer_weights.clear();
    ifstream wfile(fname);
    string kmer;
    string tmp;
    float weight;
    float threshold;
    string line;
    int index = 0; 
    while (getline(wfile, line)) {
        std::istringstream ss(line);

        if(index == 0 && (line[0]!='*')){
            cout << "Wrong file format. First line must contain threshold information" << endl;
            exit(0);
        }

        if(line[0] == '*'){
            ss >> tmp >> tmp >> threshold; 
        }else if(line[0] == '#'){
            continue;
        }else{
            std::istringstream ss(line);
            ss >> kmer >> weight;
            kmer_weights[kmer] = weight;
        }
        index++;
    }
    wfile.close();
    return threshold; 
}


// Adds up all kmer weights in the seq. Then divide by the number of kmers. 
float score_seq(const string& seq, const unordered_map<string, float>& kmer_weights, int k) { 
    float sc = 0;
    string subseq;
    int n = 0;
    for (unsigned int i = 0; i < seq.length() - k + 1; i++) {
        subseq = seq.substr(i, k);
        auto iter = kmer_weights.find(subseq);
        if (iter != kmer_weights.end()) { // kmer found
            sc += (*iter).second;
        } else { // not found
            iter = kmer_weights.find(revcomp(subseq));
            if (iter != kmer_weights.end()) { // revcomp(kmer) found
                sc += (*iter).second;
            } else { // Contains unidenfied character (i.e. not ATGC)
                cout << "Unrecognized kmer found: " << subseq << endl;
                cout<< "debugging required." << endl;
                exit(EXIT_FAILURE);
                return 0;
            }
        }
	n++;
    }
    if(n==0){
        cout<<"Length of input sequence must be larger than k"<<endl;
        exit(1);
    }

    return sc / (seq.length() - k + 1);
}
vector<float> score_seq_multigkm(const string& seq, const unordered_map<string, vector<float>>& kmer_weight_matrix, int k) {
    string subseq;
    // # gkm models
    int m = ((kmer_weight_matrix.begin())->second).size();

    vector<float> sc_v(m);
    int n = 0;
    for (unsigned int i = 0; i < seq.length() - k + 1; i++) {
        subseq = seq.substr(i, k);
        auto iter = kmer_weight_matrix.find(subseq);
        if (iter == kmer_weight_matrix.end()) { // if not  found
            iter = kmer_weight_matrix.find(revcomp(subseq));
            if (iter == kmer_weight_matrix.end()) { // even revcomp not found 
                cout << "Unrecognized kmer found: " << subseq << endl;
                cout<< "debugging required." << endl;
                exit(EXIT_FAILURE);
            }
        }
	for(int j = 0; j<m; j++){
		sc_v[j] += (iter->second)[j];
	}
        n++;
    }
    if(n==0){
        cout<<"Length of input sequence must be larger than k"<<endl;
        exit(1);
    }
    for(auto& sc : sc_v){
	sc = sc / (seq.length() - k + 1);
    }
    return sc_v;
}


// normalize gkm-SVM score to vary between 0 and 1
float posterior_transform(float sc, float mu_pos, float var_pos, float mu_neg, float var_neg){
    float sc_p;
        float P_pos = exp(-0.5*(pow(sc-mu_pos, 2))/var_pos + log(0.5) - 0.5*log(2*M_PI*var_pos));
        float P_neg = exp(-0.5*(pow(sc-mu_neg, 2))/var_neg + log(0.5) - 0.5*log(2*M_PI*var_neg));


	// in case of moderately overlapping bimodal distributions, for extremely small or large sc, both p-pos and p-neg can be close to zero. 
	// For well-seperated bimodal distributions, intermediate sc can also give p-pos and p-neg close to zero. 
	// In such cases, return 0/0.05/1
	if(P_pos < 0.001 && P_neg < 0.001){
		if(sc > mu_pos) {return 1;}
		else if(sc < mu_neg) {return 0;}
		else {return 0.5;}
	} else {
        	sc_p = (P_pos / (P_pos + P_neg));

		// keep the function as monotonic as possible without distorting it too much. 
		if(sc > mu_pos && sc_p < 0.5){
			return 1;
		}else if (sc < mu_neg && sc_p > 0.5) {
			return 0;
		}else{
			return sc_p;
		}
	}
}

// Given a sequence S of size L
// output vector v[i] stores prediction score for a window of size pred_width centered around (S[(i*slide_width + window_width) / 2])
vector<float> score_sliding_windows(string seq, int window_width, int slide_width, 
		                    tuple<int, float, float, float, float, unsigned int> post_stat,
				    const unordered_map<string, float>& kmer_weights){

    size_t dim = (seq.size() - window_width)/slide_width + 1;
    vector<float> score_vect(dim);

    int pred_width = get<0>(post_stat);
    float mu_pos = get<1>(post_stat); float var_pos = get<2>(post_stat);
    float mu_neg = get<3>(post_stat); float var_neg = get<4>(post_stat);
    unsigned kmer_size = get<5>(post_stat);

    float sc; float sc_p;
    string subseq;
    for(size_t i = 0; i < dim; i++){
    	int xo = static_cast<int>(i*slide_width - (pred_width - window_width)/2);
   	if(xo >= 0){
    		subseq = seq.substr(xo,  pred_width);
   	}else{
    		subseq = seq.substr(0,  pred_width + xo);
    	}

 	sc = score_seq(subseq, kmer_weights, kmer_size);
 	sc_p = posterior_transform(sc, mu_pos, var_pos, mu_neg, var_neg);
    	score_vect[i] = sc_p;
    }
    return score_vect;
}

// same for bigwig avg scores 
// return (genomic_windows x #experiments) matrix
vector<vector<float>> avgbw_sliding_windows(string seq, int window_width, int slide_width,
                                     vector<string> bwf_list, string chrom, int locus_start, string chain_ID){

    size_t dim = (seq.size() - window_width)/slide_width + 1;
    int N_exp = bwf_list.size();
    vector<vector<float>> score_matrix(dim, vector<float>(N_exp, 0));

    string command; 
    for(int j = 0; j<N_exp; j++){

	string bwfname = "/mnt/data0/joh27/projects/alignment_enhancer_conservation/methods/gkm_dynam/versions/whole_genome/V20/test/HBB/bigwig_outputs/" + chain_ID +"_" + chrom + "_" + to_string(locus_start) + ".bed.out" + "_" + to_string(j);
	check_file(bwfname);
        ifstream avgbwfile(bwfname);
        unsigned int index; string trash; float val; string line; 
        while (getline(avgbwfile, line)) {
           std::istringstream ss(line);
	   ss >> index >> trash >> trash >> trash >> val;
	   if(index < dim){
	   score_matrix[index][j] = val;
	   }
    	}
    }

    return score_matrix;
}




vector<vector<float>> multi_score_sliding_windows(string seq, int window_width, int slide_width,
                                    vector<tuple<int, float, float, float, float, unsigned int>> post_stat_list,
                                    const unordered_map<string, vector<float>>& kmer_weight_matrix){

    size_t dim = (seq.size() - window_width)/slide_width + 1;
    vector<vector<float>> score_matrix(dim);

    unsigned int N_gkm = ((kmer_weight_matrix.begin())->second).size(); 
    vector<float> score_vect(N_gkm);

    int pred_width = get<0>(post_stat_list[0]);
    unsigned int kmer_size = get<5>(post_stat_list[0]);
    float mu_pos; float var_pos; float mu_neg; float var_neg; 

    string subseq;
    for(size_t i = 0; i < dim; i++){
        int xo = static_cast<int>(i*slide_width - (pred_width - window_width)/2);
        if(xo >= 0){
                subseq = seq.substr(xo,  pred_width);
        }else{
                subseq = seq.substr(0,  pred_width + xo);
        }
	score_vect = score_seq_multigkm(subseq, kmer_weight_matrix, kmer_size);
        for(size_t j = 0; j < N_gkm; j++){
            auto& post_stat = post_stat_list[j];
	    mu_pos = get<1>(post_stat); var_pos = get<2>(post_stat);
            mu_neg = get<3>(post_stat); var_neg = get<4>(post_stat);
	    score_vect[j] = posterior_transform(score_vect[j], mu_pos, var_pos, mu_neg, var_neg);
        }
        score_matrix[i] = score_vect;
    }
    return score_matrix;
}


float score_bp(const string& seq, int pos, const unordered_map<string, float>& kmer_weights, int k) {
    int c1 = max(0, pos - k + 1);
    int c2 = min(static_cast<int>(seq.length()) - 1, pos + k - 1);
    string sub = seq.substr(c1, c2 - c1 + 1);
    return score_seq(sub, kmer_weights, k);
}



// used for masking 
vector<float> gen_sc_v(const string& seq, const unordered_map<string, float>& kmer_weights) {
    vector<float> bp_scores;
    int k = ((*(kmer_weights.begin())).first).length();

    for (unsigned int pos = 0; pos < seq.length(); pos++) {
        bp_scores.push_back(score_bp(seq, pos, kmer_weights, k));
    }
    return bp_scores;
}



// mask using non-ATGC ASCII
// return a boolean vector with 1 = masked 0 = not-masked
pair<string, vector<bool>> mask_seq(string seq, unordered_map<string, float>& kmer_weights, float threshold) {

    string masked_seq = seq; 
    vector<float> bp_scores = gen_sc_v(seq, kmer_weights);
    vector<bool> masked_bp(seq.length(), false);
    for (unsigned int i = 0; i < seq.length(); i++) {
        if (bp_scores[i] > threshold) {
            char c= (char) (33 + rand()%32); // special character and numbers
            masked_seq[i] = c;
	    masked_bp[i] = true;
        }
    }
    return make_pair(masked_seq, masked_bp);
}

// For seq of length N, take boolean vector of size N, then convert to percentage of size M, where M is the number of sliding windows, and ith element contain percent NOT masked of ith sliding window.
// Assume that N-w is divisible by s. Else return error. 
vector<float> bool_masked_2_percent_nmasked(vector<bool>& bool_masked, int w, int s){
	unsigned int N = bool_masked.size();
	if((N-w) % s != 0){
		cout << "debug required: in bool_masked_2_percent_nmasked, N-w must be divisible by s" << endl;
		exit(1);
	} else{

		vector<float> percent_nmasked((N-w)/s + 1, 0);
		for(unsigned int i = 0; i < percent_nmasked.size(); i++){
			int n_masked = 0; 
			for(int j = 0; j < w; j++){
				n_masked += bool_masked[i*s + j];
			}
			percent_nmasked[i] = (float)(w-n_masked) / (float)w ;
		}
		return percent_nmasked;

	}
	

}


float train_masker(const unordered_map<string, float>& kmer_weights, string rfile_name, string dir, float fraction_mask){
    check_file(rfile_name);
    ifstream rfile(rfile_name);

    int k = ((*(kmer_weights.begin())).first).length();
    vector<float> vals;  // training val

    int index;
    string line;
    while (getline(rfile, line)) { // looping through bed lines
        cout << line << endl;
        string build; string chrom;
        string start; string end; // string input to include inputs like "." for all
        string fafile_name;
        istringstream ss(line);
        ss >> build >> chrom >> start >> end;
        fafile_name = dir + build + "/" + chrom + ".fa";
        string train_seq = read_fa(fafile_name);
        transform(train_seq.begin(), train_seq.end(), train_seq.begin(), ::toupper);

        if(start == "."){
            start = "0";
        }
        if(end == "."){
            end = to_string(train_seq.size() - 1);
        }

        string subseq; int c1; int c2;
        for(int i = 0; i< (int)train_seq.size(); i++){ // generate training vals
            if(i>=stoi(start) && i<=stoi(end)){
                c1 = i - k + 1;
                c2 = i + k -1;
                if(c1>=0 && c2 < (static_cast<int>(train_seq.length())) - 1){ // full length
                    subseq = train_seq.substr(c1, c2 - c1 + 1);
                    bool unrecog_char = false;
                    for(auto elem : subseq){
                        if(!((elem == 'A') || (elem == 'C') || (elem == 'G') || (elem == 'T'))){
                            unrecog_char = true;
                            break;
                        }
                    }
                    if(!unrecog_char){ //only ATGCatgc
                        vals.push_back(score_bp(train_seq, i, kmer_weights, k));
                    }
                }
            }
        }
        cout << "#sample obtained: " <<  vals.size() << endl;

    }
    sort(vals.begin(), vals.end(), std::greater<float>());
    index = floor(fraction_mask * vals.size());

    float threshold = (vals[index] + vals[index]) / 2;
    return threshold;
}



float pcorr(vector<float> X, vector<float> Y){
   int n = X.size();
   float sum_X = 0, sum_Y = 0, sum_XY = 0;
   float squareSum_X = 0, squareSum_Y = 0;
   for (int i = 0; i < n; i++){
      sum_X = sum_X + X[i];
      sum_Y = sum_Y + Y[i];
      sum_XY = sum_XY + X[i] * Y[i];
      squareSum_X = squareSum_X + X[i] * X[i];
      squareSum_Y = squareSum_Y + Y[i] * Y[i];
   }

   float nume = (float)(n * sum_XY - sum_X * sum_Y);
   float denom = sqrt((n * squareSum_X - sum_X * sum_X) * (n * squareSum_Y - sum_Y * sum_Y)) + 0.000001;
   float corr = nume / denom;
   return corr;
}


vector<float> invert_vector(vector<float> v){
	vector<float> out(v.size(), 0);
	for(unsigned int i = 0; i < v.size(); i++){
		out[v.size() - 1 - i] = v[i];
	}
	return out;
}


vector<float> elem_prod_vectors(vector<float> v1, vector<float> v2){
	vector<float> prod(v1.size());
	for(unsigned int i; i<prod.size(); i++){
		prod[i] = v1[i]*v2[i];
	}
	return prod;
}
