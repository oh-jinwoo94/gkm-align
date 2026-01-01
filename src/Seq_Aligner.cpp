#include "header.h"
// note: matrix G and F are constructed so that main colinear conservation signals are located along the main diagonal. i.e. in case of inversion, seq2 is inverted so that seq1 and seq2 are in same orientation.
// For outputing alingment coordinate, "rel_strand_loc" information is used. 

// constructor1: without external signal (e.g. gkmSVM)
Seq_Aligner::Seq_Aligner(Matrix* G, string rel_strand_loc, string ID_2align)
:
  km_dim_1(G->dims()[0]),
  km_dim_2(G->dims()[1]),
  dm_dim_1(km_dim_1 + 1),
  dm_dim_2(km_dim_2 + 1),
  G(G),
  K(km_dim_1, km_dim_2),
  DM(dm_dim_1, dm_dim_2, numeric_limits<float>::infinity()),
  rel_strand_loc(rel_strand_loc),
  type(1),
  ID_2align(ID_2align)

{

    DM(0,0) = 0;

    // z-transform

    mean = G->compute_mean();
    var = G->compute_var();
    for(int i = 0; i < km_dim_1; i++){
        for(int j = 0; j < km_dim_2; j++){
	    K(i,j) = ((*G)(i,j) - mean) / (sqrt(var) + 0.0001);
	}
    }


}

// constructor2: with external signal
Seq_Aligner::Seq_Aligner(Matrix* G, vector<float>* ann_vector_1, vector<float>* ann_vector_2, float wF, string rel_strand_loc, string ID_2align)
:
  km_dim_1(G->dims()[0]),
  km_dim_2(G->dims()[1]),
  dm_dim_1(km_dim_1 + 1),
  dm_dim_2(km_dim_2 + 1),
  G(G),
  aV1(ann_vector_1),
  aV2(ann_vector_2),	
  K(km_dim_1, km_dim_2),
  DM(dm_dim_1, dm_dim_2, numeric_limits<float>::infinity()),
  rel_strand_loc(rel_strand_loc),
  type(2),
  ID_2align(ID_2align)

{

    DM(0,0) = 0;

    float buffer = 0.01; 
    float ann_prod;
    for(int i = 0; i < km_dim_1; i++){
        for(int j = 0; j < km_dim_2; j++){
	    ann_prod = (*aV1)[i] * (*aV2)[j];
            K(i,j) = pow((*G)(i,j), 1-wF) * (pow(ann_prod, wF) + buffer) / (1 + buffer);
        }
    }

    // z-transform
    mean = K.compute_mean();
    var = K.compute_var();
    for(int i = 0; i < km_dim_1; i++){
        for(int j = 0; j < km_dim_2; j++){
            K(i,j) = (K(i,j) - mean) / (sqrt(var) + 0.0001);
        }
    }


}




vector<tuple<int, int>> Seq_Aligner::gen_dots(float indel){
    update_DM(indel);
    vector<tuple<int, int>> dots = backtrack(indel);
    return dots;
}


void Seq_Aligner::update_DM(float indel) {

    DM(0,0) = 0;
    for (int i = 1; i < dm_dim_1; i++) {
	DM(i,0) = i * (-indel); 
    }
    for (int j = 1; j < dm_dim_2; j++) {
	DM(0,j) = j * (-indel); 
    }
    // fill in dynamic matrix
    float cand1;
    float cand2;
    float cand3;
    float gkm_sim;
    std::vector<float> v;
    float max_elem;
    for (int i = 1; i < dm_dim_1; i++) {
        for (int j = 1; j < dm_dim_2; j++) {
            gkm_sim = K(i - 1, j - 1);
           cand1 = DM(i - 1, j - 1) + gkm_sim; 
            cand2 = DM(i - 1, j) - indel;
            cand3 = DM(i, j - 1) - indel;
            v = {cand1, cand2, cand3};
            max_elem = *(max_element(begin(v), end(v)));
	    DM(i,j) = max_elem;
        }
    }
}


// alignment using K_norm
vector<tuple<int, int>>
Seq_Aligner::backtrack(float indel) {

    // returns extension coordinates (match or sub) and its matrix values
    // takes both same_strand diff_strand matrix. If diff_strand matrix, coordinate will
    // be reverse mapped.
    // i.e. orientation of the corodinates will now be consistent with the genomic orientation. 

    // allowed numerical deviation
    float epsilon = 0.01;
    if ((rel_strand_loc != "same_strand") &&
        (rel_strand_loc != "diff_strand")) {
        cout << "Unrecognized type :" + rel_strand_loc << endl;
        exit(1);
    }


    vector<tuple<int, int>> dots;
    int i = dm_dim_1 - 1;
    int j = dm_dim_2 - 1;

    tuple<int, int> dot;

    while (i > 0 || j > 0) {
        // Boolean values that encode which 'directions' are possible
        bool match = ((i > 0 && j > 0) &&
                      (fabs((DM(i - 1, j - 1) +
                             K(i - 1, j - 1)) -
                            DM(i, j)) < epsilon));
        bool indel1 = ((i > 0) &&
                       (fabs((DM(i - 1, j) - indel) - DM(i,j)) < epsilon));
        bool indel2 = ((j > 0) &&
                       (fabs((DM(i,j - 1) - indel) - DM(i,j)) < epsilon));

        if (match) {
            if (rel_strand_loc == "same_strand") {
                dot = make_tuple(i - 1, j - 1);
            } else if (rel_strand_loc == "diff_strand") {
                dot = make_tuple(i - 1, km_dim_2 - 1 - (j - 1));
            }
            dots.insert(dots.begin(), dot);
            --i;
            --j;
            num_match++;

        } else if (indel1) {
            --i;
            num_skip++;

        } else if (indel2) {
            --j;
            num_skip++;

        } else { // for debug.
            cout << "Unexpected event (debug needed)" << endl;
	    cout << ID_2align << endl;
	    (*G).save_matrix_tsv(ID_2align + ".debug.matrixG.tsv");
            K.save_matrix_tsv(ID_2align + ".debug.matrixK.tsv");
            DM.save_matrix_tsv(ID_2align + ".debug.matrixDM.tsv");

            (*G).save_matrix(ID_2align + ".debug.matrixG");
            K.save_matrix(ID_2align + ".debug.matrixK");
            DM.save_matrix(ID_2align + ".debug.matrixDM");


	    cout << "matrix mean: " << mean <<endl;
	    cout << "matrix var: " << var << endl;
            cout << DM(i,j) << endl;
	     
            cout << i << '\t' << j << endl;
            // cout<<DM[i-1][j-1]<<endl;
            cout << DM(i,j - 1) << endl;
            cout << DM(i,j - 1) - indel << endl;
            cout << fabs((DM(i,j - 1) - indel) - DM(i,j)) << endl;
            cout << (fabs((DM(i,j - 1) - indel) - DM(i,j)) < epsilon) << endl;
            exit(EXIT_FAILURE);
        }
    }
    return dots;
}



// Map DM index to genomic cooridnate (center of the window)
// coordinate = {coord1, coord2, Gval, aV1, aV2}
vector<tuple<int, int, string,  string, string>>
Seq_Aligner::dots_to_coords(
    vector<tuple<int, int>> &dots,
    int build_1_coord_start, int build_2_coord_start,
    int slide_width, int window_width, string rel_strand_loc) {
    
	
    vector<tuple<int, int, string, string,string>> coords;
    tuple<int, int, string, string, string> coord;
    int i; int j; string Gval; string aV_elem_1 = "."; string aV_elem_2 = ".";
    for (auto &dot : dots) {
        i = get<0>(dot); j = get<1>(dot);
        if(rel_strand_loc == "same_strand"){
	    if(type == 1){
		Gval = to_string((*G)(i,j));
	    } else {
	        Gval = to_string((*G)(i,j));
		aV_elem_1 = to_string((*aV1)[i]);
		aV_elem_2 = to_string((*aV2)[j]);
	    }
	}else{ // Since G and ann vectors were computed by revcomping seq_2, in case of "diff_strand", their values need to be fetched through reversed coordinate. 
            if(type == 1){
                Gval = to_string((*G)(i,km_dim_2 - j - 1));
            } else {
                Gval = to_string((*G)(i,km_dim_2 - j - 1));
		aV_elem_1 = to_string((*aV1)[i]);
                aV_elem_2 = to_string((*aV2)[km_dim_2 - j - 1]);
            }

        }
        coord = make_tuple(build_1_coord_start + (i * slide_width) +
                               window_width / 2,
                           build_2_coord_start + (j * slide_width) +
                               window_width / 2,
                           Gval, aV_elem_1, aV_elem_2);
        coords.push_back(coord);
    }
    return coords;
}



