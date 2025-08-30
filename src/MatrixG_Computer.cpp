#include "header.h"

MatrixG_Computer::MatrixG_Computer(int lmer_length,
                                 int nongap_length,
                                 int slide_width,
                                 int window_width, string seq1,
                                 string seq2, string rel_strand_loc, int threads):

    l(lmer_length),
    k(nongap_length),
    slide(slide_width),
    window(window_width),
    n_thread(threads),
    seq1(seq1),
    seq2(seq2),
    km_dim_1(1 + (seq1.length() - window_width) / slide_width),
    km_dim_2(1 + (seq2.length() - window_width) / slide_width),
    K(km_dim_1, km_dim_2, 0),
    rel_strand_loc(rel_strand_loc)
{

    // handle revcomp 
    if (rel_strand_loc == "same_strand") {
        // do nothing
    } else if (rel_strand_loc == "diff_strand") {
       seq2 = revcomp(seq2);
    } else {
        cout << "Unrecognized relative strand location: " + rel_strand_loc << endl;
        exit(1);
    }


    // SIMD vectorization: obtain all subseqs of size l
    // v[i] = seq[i:i+l] + pads
    #ifdef __AVX2__
        string pad1(32-l, '>'); // 32 for AVX2
        string pad2(32-l, '<');
    #elif defined(__SSE2__)
        string pad1(16-l, '>'); // 16 for SSE2
        string pad2(16-l, '<');
    #elif defined(__NOSIMD__)
        // No SIMD support - this code path should not be reached for -t 1
        // as it's blocked in argument parsing
        string pad1(0, '>');
        string pad2(0, '<');
    #endif
    for(unsigned int i = 0; i<seq1.length() - l + 1; i++){
	    seq1_sv_a.push_back(seq1.substr(i, l) +  pad1);
	    seq1_sv_b.push_back(seq1.substr(i, l) +  pad2);

	    seq1_pv_a.push_back( &((seq1_sv_a[i])[0]));
            seq1_pv_b.push_back( &((seq1_sv_b[i])[0]));
    }
    for(unsigned int i = 0; i<seq2.length() - l + 1; i++){
	    seq2_sv_a.push_back(seq2.substr(i, l) +  pad1);
            seq2_sv_b.push_back(seq2.substr(i, l) +  pad2);

	    seq2_pv_a.push_back( &((seq2_sv_a[i])[0]));
            seq2_pv_b.push_back( &((seq2_sv_b[i])[0]));
    }



    // initializing 3-dimenasional array
    pm_dim_1 = km_dim_1 - 1 + window / slide;
    pm_dim_2 = km_dim_2 - 1 + window / slide;

    // pre-computing values that will be repeatedly used
    for (int m = 0; m < l + 1; m++) {
        if ((l - m) >= k) {
            shared_gkm.push_back(nchoosek(l - m, k));
        } else {
            shared_gkm.push_back(0);
        }
    }


    init_norm(v1_norm, seq1_pv_a, seq1_pv_b, km_dim_1);
    init_norm(v2_norm, seq2_pv_a, seq2_pv_b, km_dim_2);
} // end of constructor 1




// computes list of vector norms for the denominator of
// <s1,s2>/sqrt(<s1,s1>)*sqrt(<s2,s2>)
void MatrixG_Computer::init_norm(vector<float> &norm, vector<char*> &seq_v_a, vector<char*> &seq_v_b, 
                                int km_dim) {
    float val;
    for (int i = 0; i < km_dim; i++) {
            val = subseqs_matrix(seq_v_a, seq_v_b, i*slide, i*slide + window - 1, i*slide, i*slide + window - 1);
            norm.push_back(val);
    }

    for (auto it = norm.begin(); it != norm.end(); it++) {
        *it = sqrt(*it);
    }
} // end of method




// given subseq indices, compute <sub1, sub2>. To be uesd in compute_matrix
float MatrixG_Computer::subseqs_matrix(vector<char*> &seq1, vector<char*> &seq2,
                                      int start_i, int end_i,
                                      int start_j,
                                      int end_j) {
    int tot = 0;
    int match;
    
    #ifdef __AVX2__
        __m256i s1, s2, ceq;
        for (int i = start_i; i <= end_i - l + 1; i++) {
            for (int j = start_j; j <= end_j - l + 1; j++) {
                s1 =  _mm256_loadu_si256((__m256i*)(seq1[i]));
                s2 =  _mm256_loadu_si256((__m256i*)(seq2[j]));
                ceq = _mm256_cmpeq_epi8(s1, s2);
                match = __builtin_popcount(_mm256_movemask_epi8(ceq));
                tot += shared_gkm[l - match];
            }
        }
    #elif defined(__SSE2__)
        __m128i s1, s2, ceq;
        for (int i = start_i; i <= end_i - l + 1; i++) {
            for (int j = start_j; j <= end_j - l + 1; j++) {
                s1 =  _mm_loadu_si128((__m128i*)(seq1[i]));
                s2 =  _mm_loadu_si128((__m128i*)(seq2[j]));
                ceq = _mm_cmpeq_epi8(s1, s2);
                match = __builtin_popcount(_mm_movemask_epi8(ceq));
                tot += shared_gkm[l - match];
            }
        }
    #elif defined(__NOSIMD__)
        // No SIMD support - this code path should not be reached for -t 1
        // as it's blocked in argument parsing
        cout << "Error: SIMD support required for sequence alignment" << endl;
        exit(1);
    #endif
 
    return static_cast<float>(tot);
} // end of method




// multithread computation
Matrix& MatrixG_Computer::compute_full_matrix(){
    vector<thread> workers;
    vector<vector<int>> allocation;
    int curr = 0;

    // leave one out for the current node. 
    for(int i = 0; i < n_thread - 1; i++){
        vector<int> vec;
        int init = curr;
        for(int j = 0; j < km_dim_1/n_thread; j++){
            vec.push_back(init + j);
            curr++;
        }
        allocation.push_back(vec);
    }
    for(int i = 0; i < n_thread-1; i++){
        workers.push_back(thread(&MatrixG_Computer::compute_matrix_rows, this, allocation[i]));
    }
    vector<int> rows;
    for(int i = curr; i < km_dim_1; i++){
        rows.push_back(i);
    }

    compute_matrix_rows(rows);
    for(auto& worker : workers){
        worker.join();
    }
    return K;

}


// compute matrix by row. Useful for multithreading
void MatrixG_Computer::compute_matrix_rows(vector<int> rows){

    if(slide>=l){ // faster version. But it uses lots of memory, so we need further memory optimization by freeing/allocating memory wisely
        int*** piece_container = new int **[pm_dim_1];
        for(unsigned int i = 0; i<rows.size(); i++){
            int row = rows[i];
            //--------------- memory handling------------------ ////
            if(i == 0){ // memory handling: initialize and allocate PM rows
                for(int a = 0; a < window/slide; a++) {
                    piece_container[row + a] = new int *[pm_dim_2];
                    for(int b = 0; b<pm_dim_2; b++) {
                        piece_container[row + a][b] = new int[4];
                        for(int c = 0; c<4; c++){
                            piece_container[row + a][b][c] = -1;
                        }
                    }
                }

            } else{ // memory handling: allocate next row in piece-container and free the row that will no longer be used
                piece_container[row + window/slide - 1] = new int *[pm_dim_2];
                for(int b = 0; b<pm_dim_2; b++) {
                    piece_container[row + window/slide - 1][b] = new int[4];
                    for(int c = 0; c<4; c++){
                        piece_container[row + window/slide - 1][b][c] = -1;
                    }
                }

                for(int j=0; j<pm_dim_2; j++){
                    delete[] piece_container[row - 1][j];
                }
                delete[] piece_container[row - 1];
            }
            // -------------------actual computation----------------//
            
            for(int j = 0; j<km_dim_2; j++){ 
                K(row, j) = compute_sliding_matrix(row,j, piece_container);
            }
            // ------------------------------------------------------//

            //memory handling: if last row, free all the remaining elements
            if(i==rows.size()-1){
                for(int r = row; r<row+window/slide; r++){
                    for(int c=0; c<pm_dim_2; c++){
                        delete[] piece_container[r][c];
                    }
                    delete[] piece_container[r];
                }
            }
        }
        delete[] piece_container;
    } else {
        for(auto row : rows) {
            float kern;
            for(int j = 0; j<km_dim_2; j++){
                kern = subseqs_matrix(seq1_pv_a, seq2_pv_b, slide*row, slide*row + window-1, slide*j, slide*j+window-1);
		K(row,j) = kern / (v1_norm[row] * v2_norm[j]);  
            }
        }
    }
}


// elements for computing matrix by sliding. 
int MatrixG_Computer::compute_piece(int pm_coord_1, int pm_coord_2, string piece_type) {
    if(piece_type == "center"){
        return subseqs_matrix(seq1_pv_a, seq2_pv_b, pm_coord_1 * slide,
                          (pm_coord_1 + 1) * slide - 1, pm_coord_2 * slide,
                          (pm_coord_2 + 1) * slide - 1);
    } else if(piece_type == "right") {
        return subseqs_matrix(seq1_pv_a, seq2_pv_b, pm_coord_1 * slide, (pm_coord_1 + 1) * slide - 1,
        (pm_coord_2 + 1) * slide - l + 1, (pm_coord_2 + 1) * slide + l - 2);

    } else if(piece_type == "down") {
        return subseqs_matrix(seq1_pv_a, seq2_pv_b, (pm_coord_1 + 1) * slide - l + 1,
                          (pm_coord_1 + 1) * slide + l - 2, pm_coord_2 * slide,
                          (pm_coord_2 + 1) * slide - 1);
    } else if (piece_type == "both") {
        return subseqs_matrix(seq1_pv_a, seq2_pv_b,(pm_coord_1 + 1) * slide - l + 1,
                                      (pm_coord_1 + 1) * slide + l - 2,
                                      (pm_coord_2 + 1) * slide - l + 1,
                                      (pm_coord_2 + 1) * slide + l - 2);
    } else {
        cout << "Debug required: unrecognized piece type: " << piece_type <<endl;
        return 0;
    }
} // end of method




// no redundant computation in sliding.
// compute row vector of the matrix matrix whose index in specified by curr_row variable. 
float MatrixG_Computer::compute_sliding_matrix(int row, int col, int*** piece_container){

    float kern = 0; 
    // now start row matrix computationi
    int pm_coord_1;
    int pm_coord_2;
    for (int a = 0; a < window / slide; a++) {
        for (int b = 0; b < window / slide; b++) {
            pm_coord_1 = row + a;
            pm_coord_2 = col + b;
            // pieces for each piece matrix elements. Use reference to update
            // the
            // target, not a copy.
            int &center = piece_container[pm_coord_1][pm_coord_2][0];
            int &right = piece_container[pm_coord_1][pm_coord_2][1];
            int &down = piece_container[pm_coord_1][pm_coord_2][2];
            int &both = piece_container[pm_coord_1][pm_coord_2][3];
            if (b == (window / slide - 1) && a < window / slide - 1) { // last column but not last row. box
                                          // with down arrow
                if (center < 0) {
                    center = compute_piece(pm_coord_1, pm_coord_2, "center");
                }
                kern += center;

                if (down < 0) {
                    down = compute_piece(pm_coord_1, pm_coord_2, "down");
                }
                kern += down;

            } else if (a == (window / slide - 1) &&
                           b < window / slide -
                                   1) { // last row but not last column. box
                                    // with right arrow
                if (center < 0) {
                    center = compute_piece(pm_coord_1, pm_coord_2, "center");
                }
                kern += center;

                if (right < 0) {
                    right = compute_piece(pm_coord_1, pm_coord_2, "right");
                }
                kern += right;

            } else if (a == (window / slide - 1) &&
                           b == (window / slide -
                                     1)) { // last row and column. box with no arrow
                if (center < 0) {
                    center = compute_piece(pm_coord_1, pm_coord_2, "center");
                }
                kern += center;

            } else { // box with both right and down arrow
                if (center < 0) {
                    center = compute_piece(pm_coord_1, pm_coord_2, "center");
                }
                kern += center;

                if (down < 0) {
                    down = compute_piece(pm_coord_1, pm_coord_2, "down");
                }
                kern += down;

                if (right < 0) {
                    right = compute_piece(pm_coord_1, pm_coord_2, "right");
                }
                kern += right;

                if (both < 0) {
                    both = compute_piece(pm_coord_1, pm_coord_2, "both"); 
                }
                kern += both;
            }
        }
    }
    return kern / (v1_norm[row] * v2_norm[col]);
} // end of method



vector<int> MatrixG_Computer::get_km_dim() {
    vector<int> dim = {km_dim_1, km_dim_2};
    return dim;
}

string MatrixG_Computer::get_rel_strand_loc() { return rel_strand_loc; }


int MatrixG_Computer::get_s() { return slide; }

int MatrixG_Computer::get_w() { return window; }

float MatrixG_Computer::sum_K() {
    float sum = 0;
    for (int i = 0; i < km_dim_1; i++) {
        for (int j = 0; j < km_dim_2; j++) {
            sum = sum + K(i,j);
        }
    }
    return sum;
}

string MatrixG_Computer::get_seq1(){
	return seq1;
}

string MatrixG_Computer::get_seq2(){
	return seq2;
}

