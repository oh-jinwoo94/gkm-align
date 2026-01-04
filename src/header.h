#ifndef HEADER
#define HEADER
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <unordered_map>
#include <string>
#include <tuple>
#include <cmath>
#include <list>
#include <limits>
#include <algorithm>
#include <unistd.h>
#include <limits.h>
#include <thread>
#include <mutex>  
#ifdef __AVX2__
#include <immintrin.h>
#elif defined(__SSE2__)
#include <emmintrin.h>
#elif defined(__ARM_NEON) 
#include <arm_neon.h>
#endif

using namespace std; 


//const vector<char> ATGC;
//const unordered_map<char, char> comp;
void help();
unsigned long long int nchoosek(int n, int k);
void check_file(string fname);


string dirname(string fname);
string basename(string fname);
string getexepath();

string revcomp(string seq);
string read_fa(string fa_fname);
string read_fa(string fa_fname, int start, int end);
pair<string, string> generate_regions_to_align(string genome_dir, string build1, string build1_chr, int build1_start, int build1_end,
                                               string build2, string build2_chr, int build2_start, int build2_end,
                                               int window, int slide);
void load_weights(unordered_map<string, float>& kmer_weights, string weight_fname);
tuple<int, float, float, float, float, unsigned int> load_post_weights(unordered_map<string, float> &kmer_weights, string fname);
tuple<int, float, float, float, float, unsigned int> load_post_weights_to_matrix(unordered_map<string, vector<float>> &kmer_weight_matrix,
                  string fname);

float load_weights_threshold(unordered_map<string, float>& kmer_weights, string weight_fname);


// assuming all in upper case. 
void filter_unknown(string& seq);

// turn to upper case, filter out unknown
string preprocess_seq(string seq);

// kmer weight based scoring 
float score_seq(const string& seq, const unordered_map<string, float>& kmer_weights, int k);
float posterior_transform(float sc, float mu_pos, float var_pos, float mu_neg, float var_neg);

// Given a sequence S of size L
// output vector v[i] stores prediction score for a window of size pred_width centered around (S[(i*slide_width + window_width) / 2])
vector<float> score_sliding_windows(string sequence, int window_width, int slide_width, tuple<int, float, float, float, float, unsigned int> post_stat, const unordered_map<string, float>& kmer_weights);

vector<vector<float>>  avgbw_sliding_windows(string seq, int window_width, int slide_width, vector<string> bwf_list, string locus_chr, int locus_start, string chain_ID); 


vector<vector<float>> multi_score_sliding_windows(string seq, int window_width, int slide_width,
                                    vector<tuple<int,  float, float, float, float, unsigned int>> post_stat_list,
                                    const unordered_map<string, vector<float>>& kmer_weight_matrix);

float score_bp(const string& seq, int pos, const unordered_map<string, float>& kmer_weights, int k);
vector<float> gen_sc_v(const string& seq, const unordered_map<string, float>& kmer_weights);
pair<string, vector<bool>> mask_seq(string seq, unordered_map<string, float>& kmer_weights, float threshold);
vector<float> bool_masked_2_percent_nmasked(vector<bool>& bool_masked, int w, int s);

float train_masker(const unordered_map<string, float>& kmer_weights, string rfile_name, string dir, float fraction_mask);


float pcorr(vector<float> X, vector<float> Y);

vector<float> invert_vector(vector<float> v);
vector<float> elem_prod_vectors(vector<float> v1, vector<float> v2);

class Matrix{
public:
    Matrix(size_t rows, size_t cols);
    Matrix(size_t rows, size_t cols, float val);
    Matrix(string ifname);
    float& operator()(size_t i, size_t j);
    float operator()(size_t i, size_t j) const;
    vector<size_t> dims();
    float compute_mean();
    float compute_var();
    int save_matrix(string fname);
    int save_matrix_tsv(string fname);
private:
    size_t n_rows;
    size_t n_cols;
    vector<float> data;
};




// input seq:
// masked nucleotides: 128 ASCII excluding ATGC (including lower-case atgc as mask)
// non-masked nucleotides: ATGC (upper-case)
class MatrixG_Computer{  // computes and stores matrix values  <seq1, seq2>

    public:

	// constructor 1: when sequence input is entered
	MatrixG_Computer(int lmer_length, int  nongap_length, int slide_width,
		int window_width, string sequence_1, string sequence_2, string rel_strand_loc, int n_thread);




	// computes list of vector norms for the denominator of <s1,s2>/sqrt(<s1,s1>)*sqrt(<s2,s2>)
	// seq_v_a and seq_v_b represent same genomic sequence but different pads for SIMD optimizaiton.
	void init_norm(vector<float> &norm, vector<char*> &seq_v_a, vector<char*> &seq_v_b, int km_dim);

        Matrix& compute_full_matrix();

	// given subseq indices, compute <sub1, sub2>. To be uesd in compute_matrix method.  
	// input is assumed to be all upper-case letters. 
        // thread_index is for multithreading (index of workers). Useful for accessing worker specific data in computing matrix
	float subseqs_matrix(vector<char*> &s1, vector<char*> &s2,
		int start_i, int end_i, int start_j, int end_j);

	int compute_piece(int pm_coord_1, int pm_coord_2, string piece_type);

        void compute_matrix_rows(vector<int> rows);

        // Used by the faster version of compute full matrix() 
        // Now takes a flat buffer and the buffer height for modulo arithmetic
        float compute_sliding_matrix(int row, int col, vector<int>& piece_buffer, int buf_height); 

        tuple<vector<float>, vector<float>> gkmsvm_predict_posterior(const unordered_map<string, float>& kmer_to_weight_1,
                                                             float mu_pos_1, float var_pos_1, float mu_neg_1, float var_neg_1,
                                                             const unordered_map<string, float>& kmer_to_weight_2,
                                                             float mu_pos_2, float var_pos_2, float mu_neg_2, float var_neg_2);

 

	vector<int> get_km_dim();

	string get_rel_strand_loc();

        void save_ksc_file(string ofname, string seq, vector<float> sc, string description);
	int get_s();
	int get_w();
	float sum_K();
        string get_seq1();
	string get_seq2();
    private: 
	int l;
	int k;
	vector<unsigned long long int> shared_gkm; // for fast gkm matrix computation
	int slide;
	int window;
        int n_thread; 
   
        string seq1;
        string seq2;


	// a and b differ in the type of pads (> or <) used. Pads are for filling up the __m128i.
	// string vectors (vector of lmers in seq)
	vector<string> seq1_sv_a; 
	vector<string> seq2_sv_a;
	vector<string> seq1_sv_b;
        vector<string> seq2_sv_b;
	// char* vector (pointing to the first letter of each lmers in seq_sv)
	vector<char*> seq1_pv_a;
        vector<char*> seq2_pv_a;
	vector<char*> seq1_pv_b;
        vector<char*> seq2_pv_b;


	int km_dim_1;
	int km_dim_2;
	int pm_dim_1;
	int pm_dim_2;

        Matrix K; // matrix matrix

 
	// vector<float> v1_norm;
	// vector<float> v2_norm;
        vector<float> v1_inv_norm; // 1/norm (saves computation time)
        vector<float> v2_inv_norm;
	string rel_strand_loc; // currently: same_strand, diff_strand. Later include RC=TRUE


}; // end of MatrixG_Computer class



// class Seq_Aligner
// 'K_norm' is the Z-normalized gkm-matrix used for dynamic programming, while K is the orignal gkm matrix matrix.
class Seq_Aligner{  

    public:
	// constructor
	Seq_Aligner(Matrix* G,  string rel_strand_loc, string ID_2align);
        Seq_Aligner(Matrix* G, vector<float>* ann_vector_1, vector<float>* ann_vector_2, float wF,  string rel_strand_loc,  string ID_2align);
	//destructor
//	~Seq_Aligner();


        vector<tuple<int, int>> gen_dots(float indel);

	void update_DM(float indel);

	vector<tuple<int,int>> backtrack(float indel);



	// Map DM index to genomic cooridnate (center of the window)
	vector<tuple<int, int,string, string, string>> dots_to_coords(vector<tuple<int, int>>& dots,
		int build_1_coord_start, int build_2_coord_start, int slide_width, int window_width, string rel_strand_loc);



    private:
	int km_dim_1;
        int km_dim_2;
	int dm_dim_1;
        int dm_dim_2;
	Matrix* G;
	vector<float>* aV1;
	vector<float>* aV2;
	Matrix K;
        Matrix DM; // dynamic matrix 
        string rel_strand_loc;
	int type;
	float mean;
	float var;
        string ID_2align;
	
};

class Mapper{

    public:
        Mapper(string coord_file, string qbuild, vector<string> qchrom_list);
        vector<string> identify_chains(string chr, int start, int end);
        vector<tuple<int, int, string>>& get_coords(string chain_id);
        tuple<int, int, string> search_query(vector<tuple<int, int, string>>& coord_list, int qstart, int qend);
        string get_tchrom(string chainID);
        string get_tbuild();
        void map_beds();
        void save_mapfile(string ofname);


    private:
        int qdim; // In the 2D coordinates, dimension of query build. (first or second)
        string qbuild;
        string tbuild;        
        // chrom -> {chain_id, range_begin, range_end}
        // used to find a list of relevant chains that match query chromosome
        map<string, vector<tuple<string, int, int>>> chrom_to_cinfo;

        // chain_ID -> vector of {coord1, coord2, gkm_sim}. coord1: query. coord2: target
        // Given a list of chains, fetch the list of coordinates
        map<string, vector<tuple<int, int, string>>> chainID_to_coords;

        // chain_ID -> target chromosome.
        // Used for outputing result
        map<string, string> chainID_to_tchrom;
      
};
#endif
