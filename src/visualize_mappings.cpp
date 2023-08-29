#include "header.h"

vector<float> track_shared_gkmer_contribution(string which_seq, string seq1, string seq2, int lmer_length, Matrix& G, bool);
void save_cons_file(string ofname, string seq, vector<float> sc, string description);
vector<float> avg_deltaSVM(string& seq, unordered_map<string, float>& kmer_weights);
vector<char> ATGC = {'A', 'T', 'G', 'C'};

int main(int argc, char *argv[]) {
    int lmer_length = 11;
    int non_gap = 7;
    bool colinear = false; // visualize colinearly conserved gapped kmers only. (i.e. TFBS in the same order in human and mouse enhancers)
    string gdir = "NULL";
    string rmfile_name = "NULL";
    string ifname;
    string model_list_fname = "NULL";
    string ofile_directory = "";
    bool output_G_matrix = false;
    string of_prefix = "NULL";    
    char opt;
    while ((opt = getopt(argc, argv, ":l:k:W:g:d:n:o:cGh")) != -1) {
        switch (opt) {
        case 'l':
            lmer_length = atoi(optarg);
            break;
        case 'k':
            non_gap = atoi(optarg);
            break;
	case 'c':
	    colinear = true;
	    break;
        case 'G':
            output_G_matrix = true;
            break;
	case 'd':
	    gdir = string(optarg);
	    break;
        case 'g':
            rmfile_name = string(optarg);
            break;

	case 'W':
	    model_list_fname = string(optarg);
	    break;
        case 'h':
            cout << "description" << endl;
            exit(0);
	case 'n':
	    of_prefix = string(optarg);
	    break; 
        case 'o':
            ofile_directory = string(optarg);
            if (*(prev(ofile_directory.end())) != '/') {
                ofile_directory = ofile_directory + '/';
            }
            break;

        case ':':
            printf("Missing argument for %c\n", optopt);
            exit(1);

        case '?':
            cout << "Unknown option: " << char(optopt) << endl;
            exit(1);
        }
    }

    if (optind + 1 != argc) {
        cout << "One positional argument is required: input file name"
                 << endl;
            exit(1);
    } else {
            ifname = argv[optind];
    }

    if (gdir == "NULL"){
	    cout << "genome directory (containing directories for both query and target genome builds)" << endl;
	    exit(1);
    }
 


    cout << "\n-- Parameters --" << endl;
    cout << "l-mer length: " << lmer_length << " (default: 11)" << endl;
    cout << "number of informative base pairs: " << non_gap
         << " (default: 7)" << endl;

    if(colinear){cout << "colinear common gapped kmer visualization" << endl;}
    else{cout << "visualizing all common gapped kmers" << endl;};


    // load weight files, if provided 
    vector<tuple<string,unordered_map<string, float>>> model_list;
    string model_id; string model_fn;
    unordered_map<string, float> kmer_weights;
    if(model_list_fname != "NULL" ){
        check_file(model_list_fname);
        ifstream ml_file(model_list_fname);
        string line;
        while (getline(ml_file, line)) {
            istringstream ss(line);
            ss >>  model_id >> model_fn;
            check_file(model_fn);
            load_weights(kmer_weights, model_fn);
            model_list.push_back(make_tuple(model_id, kmer_weights));
        }
    }

    // load maskers, if provided
    unordered_map<string, float> masker_weights_1;
    unordered_map<string, float> masker_weights_2;
    float mthreshold_1;
    float mthreshold_2;
    if(rmfile_name != "NULL"){
	ifstream rm_file (rmfile_name);
        string rmodel_fname1;
        string rmodel_fname2;
        getline(rm_file, rmodel_fname1);
        getline(rm_file, rmodel_fname2);

        mthreshold_1 = load_weights_threshold(masker_weights_1, rmodel_fname1);
        mthreshold_2 = load_weights_threshold(masker_weights_2, rmodel_fname2);
    }

    check_file(ifname);
    ifstream ifile(ifname);
    string line; 
            //mm10    chr2    152437520       152437820       -->     hg38    chr20   279005  279304    qDnase:150_mDnase:120_gkmsim_0.2_mgkmSVM_0.98_reg_152.3

    string build1; string chr1; int start1; int end1;
    string build2; string chr2; int start2; int end2;
    string cons_info;
    string tmp;
    int index = 0;
     while (getline(ifile, line)) {
         istringstream ss(line);
         ss >> build1 >> chr1 >> start1 >> end1 >>  tmp
                     >> build2 >> chr2 >> start2 >> end2 >> cons_info;

         if(chr2 == "not-mapped"){
             continue;
         }
         string build1_fname = gdir + "/" + build1 + "/" + chr1 + ".fa";
         string build2_fname = gdir + "/" + build2 + "/" + chr2 + ".fa";
         string seq_1 = read_fa(build1_fname, start1, end1);
         string seq_2 = read_fa(build2_fname, start2, end2);
         seq_1 = preprocess_seq(seq_1);
         seq_2 = preprocess_seq(seq_2);
	 string mseq_1 = seq_1;
	 string mseq_2 = seq_2;
	 if(rmfile_name != "NULL"){
                auto masking_pair_1 = mask_seq(seq_1, masker_weights_1, mthreshold_1);
                auto masking_pair_2 = mask_seq(seq_2, masker_weights_2, mthreshold_2);
                mseq_1 = masking_pair_1.first;
                mseq_2 = masking_pair_2.first;
	 }
	 // first test if the two sequences are in the same or different strands by computing seq similairty in both directions. 
	 MatrixG_Computer kc_tmp1(lmer_length, non_gap, 1,  seq_1.length(),
			 seq_1, seq_2, "same_strand", 1);
	 Matrix& G1_tmp = kc_tmp1.compute_full_matrix(); 
	
	 MatrixG_Computer kc_tmp2(lmer_length, non_gap, 1,  seq_1.length(),
                               seq_1, revcomp(seq_2), "same_strand", 1);
         Matrix& G2_tmp = kc_tmp2.compute_full_matrix();
   
	 if(G1_tmp(0,0) < G2_tmp(0,0)){
		seq_2 = revcomp(seq_2);
 	 	mseq_2 = revcomp(mseq_2);
	 }


	 // now compute matrix to track where shared gapped kmers come from
	 MatrixG_Computer kc(lmer_length, non_gap, 1, lmer_length, mseq_1, mseq_2, "same_strand", 1);
	 Matrix& G = kc.compute_full_matrix();


	 if(colinear){  
		// align within enahcner pairs, and set matrix elements that are included in the path to 0
		Seq_Aligner sa(&G, "same_strand", "");
                vector<tuple<int, int>> aligned_dots = sa.gen_dots(0);
		auto dims = G.dims();
		Matrix G_bools = Matrix(dims[0], dims[1],0);
		for(auto dot : aligned_dots){
			G_bools(get<0>(dot), get<1>(dot)) = 1;
		}

		for(unsigned int i = 0; i < dims[0]; i++){
		    for(unsigned int j = 0; j < dims[1]; j++){
	                if(G_bools(i,j) == 0){
				G(i,j) = 0;
			}
		    }
		}

	 }

	 if(output_G_matrix){
	     G.save_matrix(of_prefix + "_" +to_string(index) + ".matrixG");
	
	 
	 }

	// diff normalization is applied by whether it's colinear or not 
         vector<float> seq1_cons_sc = track_shared_gkmer_contribution("1", mseq_1, mseq_2, lmer_length, G, colinear);
         vector<float> seq2_cons_sc = track_shared_gkmer_contribution("2", mseq_1, mseq_2, lmer_length, G, colinear);

	 if(of_prefix == "NULL"){
             of_prefix = ifname;
	 }
         string ofname_q = ofile_directory +  of_prefix + "_q" + to_string(index) + ".cons_sc";
         string ofname_t = ofile_directory + of_prefix + "_t" + to_string(index) + ".cons_sc";
         save_cons_file(ofname_q, seq_1, seq1_cons_sc, line);
         save_cons_file(ofname_t, seq_2, seq2_cons_sc, line);

         string id_q = build1 + "/" +  chr1 + ":" + to_string(start1) + "-" + to_string(end1) + "/";
         string id_t = build2 + "/" +  chr2 + ":" + to_string(start2) + "-" + to_string(end2) + "/";

         string ofname = ofile_directory + of_prefix + "_" + to_string(index) + ".png";

         string command;
         // assumes python script is in the executable's directory.
         if(model_list_fname == "NULL"){
             command = "python3 " + dirname(getexepath()) + "/" +  "logo_shared_gkmer.py "
                                                   + ofname_q + " " + ofname_t + " " + id_q + " " + id_t + " " + cons_info
                                                   + " " + ofname;
         } else { // call a different python script after computing scores.
             string id;
             vector<float> bp_scores;
             vector< tuple<string, vector<float>> > score_matrix_q;
             for(auto& model : model_list){
                 id = get<0>(model);
                 bp_scores = avg_deltaSVM(seq_1, get<1>(model));
                 score_matrix_q.push_back(make_tuple(id, bp_scores));
             }

             vector< tuple<string, vector<float>> > score_matrix_t;
             for(auto& model : model_list){
                 id = get<0>(model);
                 bp_scores = avg_deltaSVM(seq_2, get<1>(model));
                 score_matrix_t.push_back(make_tuple(id, bp_scores));
             }

	     string gkmsvm_weight_ofile_q = ofile_directory +  ifname + "_q" + to_string(index) + ".avg_dSVM";
             string gkmsvm_weight_ofile_t = ofile_directory +  ifname + "_t" + to_string(index) + ".avg_dSVM";

	     // output gkmsvm score on query sequence
             ofstream gkmsvm_ofile_q(gkmsvm_weight_ofile_q);
             gkmsvm_ofile_q << "gkmsvm_models" << '\t';
             for(auto model : score_matrix_q){
                 gkmsvm_ofile_q << get<0>(model) << '\t';
             }
             gkmsvm_ofile_q << endl;

             for(unsigned int i = 0; i < seq_1.size(); i++){
                 gkmsvm_ofile_q << seq_1[i] << '\t';
                 for(auto model : score_matrix_q){
                     gkmsvm_ofile_q << get<1>(model)[i] << '\t';
                 }
                 gkmsvm_ofile_q << endl;
             }
             gkmsvm_ofile_q.close();
	     


	     ofstream gkmsvm_ofile_t(gkmsvm_weight_ofile_t);
             gkmsvm_ofile_t << "gkmsvm_models" << '\t';
             for(auto model : score_matrix_t){
                 gkmsvm_ofile_t << get<0>(model) << '\t';
             }
             gkmsvm_ofile_t << endl;

             for(unsigned int i = 0; i < seq_2.size(); i++){
                 gkmsvm_ofile_t << seq_2[i] << '\t';
                 for(auto model : score_matrix_t){
                     gkmsvm_ofile_t << get<1>(model)[i] << '\t';
                 }
                 gkmsvm_ofile_t << endl;
             }
             gkmsvm_ofile_t.close();

             command = "python3 " + dirname(getexepath()) + "/" +  "logo_shared_gkmer_and_deltaSVM.py "
                                                   + ofname_q + " " + gkmsvm_weight_ofile_q + " " + id_q + " "
                                                   + ofname_t + " " + gkmsvm_weight_ofile_t + " " + id_t + " " + cons_info
                                                   + " " + ofname;

         }
             system(command.c_str());
             cout<<"command: " << command<<endl;
             index++;

    }

    

    
    return 0;
}



vector<float> track_shared_gkmer_contribution(string which_seq, string seq1, string seq2, int l,  Matrix& G, bool colinear){
    vector<float> sums;
    vector<size_t> dims = G.dims();
    int seqL;
    int w;
    if(which_seq == "1"){
        seqL = seq1.size();
	w = seq2.size();
        for(int i = 0; i < (int)dims[0]; i++){
            float sum = 0;
            for(int j = 0; j < (int)dims[1]; j++){
                sum = sum + G(i,j);
            }
            sums.push_back(sum);
        }

    }else if(which_seq == "2"){ // col sum
        seqL = seq2.size();
	w = seq1.size();
        for(int j = 0; j < (int)dims[1]; j++){
            float sum = 0;
            for(int i = 0; i < (int)dims[0]; i++){
                sum = sum + G(i,j);
            }
            sums.push_back(sum);
        }
    }else{
        cout<<"Unrecognized input: " << which_seq << endl;
        exit(1);
    }

    vector<float> vo;
    int lower_bound; int upper_bound; float  val;
    for(int i = 0; i<=seqL-1; i++){
        lower_bound = min(max(0, i-l+1), seqL-l);
        upper_bound = max(min(seqL-l, i), 0);
        val = 0;
        for(int j = lower_bound; j<=upper_bound; j++){
            val = val + sums[j];
        }
//normalize to 1 
	if(colinear){
        	vo.push_back(val/l);
	}else{
		vo.push_back(val/(l*(w-l+1)));
	}

    }
    return vo;
}


void save_cons_file(string ofname, string seq, vector<float> sc, string description){
    ofstream ofile(ofname);
    ofile << "# " << endl;
    ofile << "# " << description << endl;
    ofile << "# " << endl;

    ofile << "pos" << '\t' << 'A' << '\t' << 'C' << '\t'<< 'G' << '\t' <<'T' << endl;

    for(unsigned int i = 0; i<seq.size(); i++){
        vector<float> orow = {0,0,0,0};
        if(seq[i] == 'A'){
            orow[0] = sc[i];
        }else if(seq[i] == 'C'){
            orow[1] = sc[i];
        }else if(seq[i] == 'G'){
            orow[2] = sc[i];
        }else if(seq[i] == 'T'){
            orow[3] = sc[i];
        }else{
            cout<<"Unexpected character: " << seq[i] <<endl;
        }

        ofile << i << '\t' << orow[0] << '\t' << orow[1] << '\t' << orow[2] << '\t' << orow[3] <<  endl;
    }
    ofile.close();
}

vector<float> avg_deltaSVM(string& seq, unordered_map<string, float>& kmer_weights) {
    transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    vector<float> bp_scores;
    int k = ((*(kmer_weights.begin())).first).length();
    char original_nuc; float init_sc; float sum;
    for (unsigned int pos = 0; pos < seq.length(); pos++) {
        original_nuc = seq[pos];
        init_sc = score_bp(seq, pos, kmer_weights, k);
        sum = 0;
        for(char SNV : ATGC){
            if(SNV != original_nuc){
                seq[pos] = SNV;
                sum = sum + (score_bp(seq, pos, kmer_weights, k) - init_sc);
            }
        }
        seq[pos] = original_nuc;
        bp_scores.push_back(sum/((ATGC.size()-1)));
    }
    return bp_scores;
}

