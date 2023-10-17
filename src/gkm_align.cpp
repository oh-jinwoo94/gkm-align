#include "header.h"
mutex mtx;

int main(int argc, char *argv[]) {
    srand(1);
    // default parameters
    int lmer_length = 11;
    int non_gap = 7;
    int slide_width = 20;
    int window_width = 300;
    float indel = 0;


    unsigned int n_thread = 1;

    string input_file_type = "NULL";
    // genome directory
    string dir = "NULL";


    // gkm-repeat masker file. contains two lines, each containing gkm-repeat masker file name for build1 and build2
    string rmfile_name = "NULL";

    string annotation_fname = "NULL";
    string annotation_type = "NULL";
    // weighted alignment
    float wF;

    // for bigwig weighting 
    string bwavg_bin_fname; // the bigwigaverageoverbed path
    vector<string> bwf_list_1;
    vector<string> bwf_list_2;
    

    string ifile_name;
    string ofile_directory = "";
    string ofname_prefix = "NULL";
    string cfile_name = "NULL";


    bool output_G_matrix = false;
    bool overwrite = false;

    // for mapping
    string qbuild = "NULL";
    bool multiple_mapping = false;
    bool unique_mapping = false;

    // post gkmsvm weight for alignment
    char opt;
    while ((opt = getopt(argc, argv, ":l:k:s:w:W:i:d:t:g:c:q:p:n:o:mubGOh")) != -1) {
        switch (opt) {
        case 'l':
            lmer_length = atoi(optarg);
            break;
        case 'k':
            non_gap = atoi(optarg);
            break;
        case 's':
            slide_width = atoi(optarg);
            break;
        case 'w':
            window_width = atoi(optarg);
            break;
        case 'i':
            indel = stof(optarg);
            break;
        case 'd':
            dir = string(optarg);
            break;
        case 'g':
            rmfile_name = string(optarg);
            break;
	case 'W': // input type: 0.1,filename
	    //stringstream ss("123,456");//string(optarg));
	    {
	    stringstream ss(optarg);
	    annotation_type = "gkm-SVM";
            string tmp_val;
            getline(ss, tmp_val, ',');
            wF = stof(tmp_val);
            getline(ss, annotation_fname, ',');

	    if(wF<0 || wF>1){
		    cout << "wF must be between 0 and 1" << endl;
		    exit(1);
	    }
	    cout << "Weighted Alignment" << endl;
            cout << "File name: " << annotation_fname << endl;
	    cout << "Weight: " << wF << endl;

	//    help();
            break;
	    }
        case 'c':
            cfile_name = string(optarg);
            break;
	case 'm':
            multiple_mapping = true;
            break;
	case 'u':
            unique_mapping = true;
            break;
        case 't':
            input_file_type = optarg;
            if (input_file_type != "1" &&
                input_file_type != "2") {
                cout << "1 (genomic intervals) 2 (chain for mapping)"
                     << endl;
		help();
                exit(1);
            }
            break;
        case 'p':
            n_thread = atoi(optarg);
            break;

        case 'q':
            qbuild = string(optarg);
            break;
        case 'o':
            ofile_directory = string(optarg);
            if (*(prev(ofile_directory.end())) != '/') {
                ofile_directory = ofile_directory + '/';
            }
            break;
        case 'n':
	    ofname_prefix = string(optarg);
	    break;

        case 'G':
            output_G_matrix = true;
            break;

        case 'O': // currently only for matrix computation.
            overwrite = true;
            break;
        case 'h':
	    help();
            exit(0);

        case ':':
            printf("Missing argument for %c\n", optopt);
	    help();
            exit(1);

        case '?':
            cout << "Unknown option: " << char(optopt) << endl;
	    help();
            exit(1);
        }
    }

    if (input_file_type == "NULL") {
        cout << "option -t (input file type) is required" << endl;
	help();
        exit(1);
    }


    if (input_file_type == "1") { // genomic intervals
        if (dir == "NULL") {
            cout << "For input type 1, option -d (genome file location) is "
                    "required"
                 << endl;
	    help();
            exit(1);
        }
        if (rmfile_name == "NULL") {
            cout << "For input type 1, option -g (genome kmer weight file for masker) is "
                    "required"
                 << endl;
	    help();
            exit(1);
        }


        cout << "\n-- Parameters --" << endl;
        cout << "insertion/deletion penalty: " << indel << endl;
        cout << "l-mer length: " << lmer_length << " (default: 19)" << endl;
        cout << "number of informative base pairs: " << non_gap
             << " (default: 7)" << endl;
        cout << "slide width: " << slide_width << " (default: 20)" << endl;
        cout << "window width: " << window_width << " (default: 80)" << endl;
        cout << "number of parallel threads: " << n_thread << " (default: 1)" << endl;

        if (output_G_matrix) {
            cout << "Matrix matrix G (general sequence similarity) will be generated for each alignment"
                 << endl;
        }


        if (optind + 1 != argc) {
            cout << "One positional argument is required: input file name"
                 << endl;
	    help();
            exit(1);
        } else {
            ifile_name = argv[optind];
        }

        check_file(ifile_name);
        ifstream ifile(ifile_name);



        // output file info
        string init_dir = dirname(ifile_name);
        if (ofile_directory == "") { // if not specified, save it where the input file is located
            ofile_directory = init_dir;
        }
	string coord_ofname;
	if(ofname_prefix != "NULL"){
            coord_ofname = ofile_directory + ofname_prefix + ".coord";
	}else{
	    coord_ofname = ofile_directory + basename(ifile_name) + ".coord"; 
	}


        // prep repeat masker
	ifstream rm_file (rmfile_name);
	string rmodel_fname1;
        string rmodel_fname2;
	getline(rm_file, rmodel_fname1);
        getline(rm_file, rmodel_fname2);

        cout << "Reading masker model: " << rmodel_fname1 <<endl;
        cout << "Reading masker model: " << rmodel_fname2 <<endl;

	
        unordered_map<string, float> masker_weights_1;
        unordered_map<string, float> masker_weights_2;

        const float mthreshold_1 = load_weights_threshold(masker_weights_1, rmodel_fname1);
	const float mthreshold_2 = load_weights_threshold(masker_weights_2, rmodel_fname2);


	tuple<int, float, float, float, float, unsigned int> post_stat_1;
	unordered_map<string, float> kmer_weights_1; // encoding reg vocab
        tuple<int, float, float, float, float, unsigned int> post_stat_2;
        unordered_map<string, float> kmer_weights_2; // encoding reg vocab


	unordered_map<string, vector<float>> kmer_weight_matrix;

	if(annotation_type == "gkm-SVM"){
	// posterior weight file for weighted alignment.
        // w    300
        // mu_pos  -0.4698601590268965
        // var pos 0.0002796103768422355
       // mu_neg  -0.5784510419876551
        // var_neg 0.0010991364946907067
        // CCCCGATATAC     -0.6253484
	    ifstream ifile(annotation_fname);
	    string model_fname;
	    getline(ifile, model_fname);
	    post_stat_1 = load_post_weights(kmer_weights_1, model_fname);
	    getline(ifile, model_fname);
	    post_stat_2 = load_post_weights(kmer_weights_2, model_fname);

        }
	

    //input format
     //some header
	//hg38    chr12   98635113        98655141        mm10    chr10   91072770        91092940        same_strand     chain1
	//
        // count the number of interval lines (excluding the header)
        // This determines how parallel processing will work
        // If nlines>p, each threads take different interval lines
        // else, multiple threads are used to computed a single gkm matrix 
        vector<tuple<string, string, int, int, string, string, int, int, string, string>> all_intervals; 
        string line;
        getline(ifile, line); // skip the header
        while(!ifile.eof()) {
	    getline(ifile, line);
            if(line.length() == 0 ){continue;} // skip empty lines
            istringstream ss(line);
      
            string build1;
            string build1_chr;
            int build1_start;
            int build1_end;

            string build2;
            string build2_chr;
            int build2_start;
            int build2_end;

            string rel_strand_loc; // whether build_1 and build_2's genetic elements are located on the same strand.
            string interval_ID; // e.g. chain id. Doesn't need to be unique

            ss >> build1 >>  build1_chr >> build1_start >>
                build1_end >> build2  >> build2_chr >> build2_start
               >> build2_end  >> rel_strand_loc >> interval_ID;

	    // if negative coordinate provided (possibily by mistake, use max(coord, 0))
	    // going over boundary isn't a problem since read_fa takes that into acccount
	    build1_start = build1_start > 0 ? build1_start : 0;
	    build2_start = build2_start > 0 ? build2_start : 0;

            all_intervals.push_back(make_tuple(build1, build1_chr, build1_start, build1_end,
                                           build2, build2_chr, build2_start, build2_end,
                                           rel_strand_loc, interval_ID));
        }

        // determine thread parameters, and partition the intervals accordingly 
        vector<vector<tuple<string, string, int, int, string, string, int, int, string, string>>> interval_partition;
        int threads_per_matrix;
        if(all_intervals.size() < n_thread){ // |partition| = 1
            threads_per_matrix = n_thread;
            interval_partition.push_back(all_intervals);
        }else{ // |partition| = n_thread
            threads_per_matrix = 1;
            for(unsigned int i = 0; i<n_thread; i++){ // insert empty vectors
                interval_partition.push_back({});
            }

            for(unsigned int i = 0; i<all_intervals.size(); i++){ // partition
                (interval_partition[i%n_thread]).push_back(all_intervals[i]);
            }
            
        }



	// lambda function for parallel threads Matrix F(dim_1, dim_2);
        // to be filled through lambda functions of the threads 
        // chain_ID -> {coordinates, chain_info}
	// coordinates = {coord1, coord2, Gval, Fval}  in case of unweighted alignment, Fval = "."
	// chain_info just contains build, chrom, rel_strand_loc
        map<string, pair<vector<tuple<int, int, string, string, string>>, vector<string>>> chainID_to_coords;
        // lambda function to distribute over multiple threads   
        auto align_intervals = [&](vector<tuple<string, string, int, int, string, string, int, int, string, string>> intervals){ 
           // input: intervals
            static int progress = 0;
            vector<tuple<int, int, string, string, string>> coords; // to be filled after alignment 
            vector<string> chain_info; // For each chain, info on build, chrom, and direction
            string chain_ID; // Intervals from the same chain have the same ID.

            for(auto it = intervals.begin(); it != intervals.end(); it++){ 
                // redefine these variables to avoid thread collisions. 
                string build1 = get<0>(*it);
                string build1_chr = get<1>(*it);
                int build1_start =  get<2>(*it);
                int build1_end = get<3>(*it);

                string build2 = get<4>(*it);
                string build2_chr = get<5>(*it);
                int build2_start = get<6>(*it);
                int build2_end =  get<7>(*it);
                string rel_strand_loc = get<8>(*it); // whether build_1 and build_2's genetic element are located on the same strand.
                string chain_ID = get<9>(*it); // e.g. chain id. Doesn't need to be unique

		if(!(rel_strand_loc =="same_strand" || rel_strand_loc == "diff_strand")){
			cout << "9'th element in the input .2align file must be either 'same_strand' or 'diff_strand'" <<endl;
			exit(1);
		}

		// unique identifier for debugging 
		string ID_2align = build1 + "-" + build1_chr + "-" + to_string(build1_start) + "-" + to_string(build1_end) + "-" +  build2 + "-" + build2_chr + "-" + to_string(build2_start) + "-" + to_string(build2_end) + "-" + rel_strand_loc; 

                chain_info = {build1, build1_chr, build2, build2_chr, rel_strand_loc}; 

                string matrix_pref = ofile_directory + build1 + "-" + build1_chr +
                                "-" + to_string(build1_start) + "-" +
                                to_string(build1_end) + "-" + build2 + "-" +
                                build2_chr + "-" + to_string(build2_start) +
                                "-" + to_string(build2_end) + "-" + rel_strand_loc;

                ifstream ifile(matrix_pref  + ".matrixG");

		// fetch seq 
		string seq_1; string seq_2;
		string mseq_1; string mseq_2; // masked sequences 
		vector<bool> bool_masked_1; vector<bool> bool_masked_2;
		vector<float> percent_nmasked_1; vector<float> percent_nmasked_2;

                pair<string, string> seqs =generate_regions_to_align(dir,
                                    build1, build1_chr, build1_start, build1_end,
                                    build2, build2_chr, build2_start, build2_end,
                                               window_width, slide_width);

                seq_1 = get<0>(seqs);
                seq_2 = get<1>(seqs);

		if((int)seq_1.size() < window_width || (int)seq_2.size() < window_width){continue;}

		// repeat masking
		auto masking_pair_1= mask_seq(seq_1, masker_weights_1, mthreshold_1);
                auto masking_pair_2 = mask_seq(seq_2, masker_weights_2, mthreshold_2);
		mseq_1 = masking_pair_1.first;
		bool_masked_1 = masking_pair_1.second;
		mseq_2 = masking_pair_2.first;
                bool_masked_2 = masking_pair_2.second;


                percent_nmasked_1 = bool_masked_2_percent_nmasked(bool_masked_1, window_width, slide_width);
                percent_nmasked_2 = bool_masked_2_percent_nmasked(bool_masked_2, window_width, slide_width);

       
		// dimension for gkmsim matrix. 
		// To ensure matrix dimension consistency, -w and -s paramters must be kept consistent for de-novo and precomputed matrices.
		size_t dim_1 = (seq_1.size() - window_width)/slide_width + 1;
		size_t dim_2 = (seq_2.size() - window_width)/slide_width + 1;

       	       
		// For weighted alignment
		vector<float> ann_vector_1;
		vector<float> ann_vector_2;
                if(annotation_type == "gkm-SVM"){
		    ann_vector_1 = score_sliding_windows(seq_1, window_width, slide_width, post_stat_1, kmer_weights_1);
		    if(rel_strand_loc == "same_strand"){
                    	ann_vector_2 = score_sliding_windows(seq_2, window_width, slide_width, post_stat_2, kmer_weights_2);
		    }else{
	                ann_vector_2 = score_sliding_windows(revcomp(seq_2), window_width, slide_width, post_stat_2, kmer_weights_2);
		    }

	            ann_vector_1 = elem_prod_vectors(ann_vector_1, percent_nmasked_1); // to correctly weigh repeats .. 
	            ann_vector_2 = elem_prod_vectors(ann_vector_2, percent_nmasked_2);
		}

                if(!overwrite && ifile.good()){ // using pre-exising matrix G
                    cout << "Relevant matrix file already exists in the ouptut directory: "<< matrix_pref  + ".matrixG" <<endl;
                    cout << "The matrix file will be used for alignment" << endl;
                    cout << "Use option -O to overwrite" << endl;
		    Matrix G = Matrix(matrix_pref  + ".matrixG");
		    auto dims = G.dims();
		    if(dims[0] != dim_1 || dims[1] != dim_2){
			    cout << "mismatch between pre-computed matrix dimension and input dimension. Check paramters" << endl;
			    exit(1);
		    } 
                    if(annotation_type == "NULL"){ // unweighted alignment
                    	Seq_Aligner sa(&G, rel_strand_loc, ID_2align);
			vector<tuple<int, int>> aligned_dots = sa.gen_dots(indel);
                        coords =  sa.dots_to_coords(aligned_dots, build1_start, build2_start,
                                  slide_width, window_width, rel_strand_loc);

		    } else { 	 // weighted alignment
			Seq_Aligner sa(&G, &ann_vector_1, &ann_vector_2, wF, rel_strand_loc, ID_2align);
			vector<tuple<int, int>> aligned_dots = sa.gen_dots(indel);
                        coords =  sa.dots_to_coords(aligned_dots, build1_start, build2_start,
                                  slide_width, window_width, rel_strand_loc);
		    }

                } else { // de novo matrix computation
                    MatrixG_Computer kc(lmer_length, non_gap, slide_width, window_width,
                               mseq_1, mseq_2, rel_strand_loc, threads_per_matrix);

		    if((dim_1 == 1) || (dim_2 ==1)){
                        cout << "Alignment skipped because gkm matrix width or height is too small"<< endl;
			mtx.lock();
                        progress++;
                        mtx.unlock();
                        continue;
                    }

                    Matrix& G = kc.compute_full_matrix();
                    if (output_G_matrix) {
			G.save_matrix(matrix_pref + ".matrixG");
                    }

                    if(annotation_type == "NULL"){ // unweighted alignment
                        Seq_Aligner sa(&G, rel_strand_loc, ID_2align);
                        vector<tuple<int, int>> aligned_dots = sa.gen_dots(indel);
                        coords =  sa.dots_to_coords(aligned_dots, build1_start, build2_start,
                                  slide_width, window_width, rel_strand_loc);

                    } else { // weighted alignment

                        Seq_Aligner sa(&G, &ann_vector_1, &ann_vector_2, wF, rel_strand_loc, ID_2align);

                        vector<tuple<int, int>> aligned_dots = sa.gen_dots(indel);

                        coords =  sa.dots_to_coords(aligned_dots, build1_start, build2_start,
                                  slide_width, window_width, rel_strand_loc);


                    }
                } // end of if-else

                mtx.lock();
                if(chainID_to_coords.find(chain_ID) == chainID_to_coords.end()){ // new chain
                    chainID_to_coords[chain_ID] = make_pair(coords, chain_info);   
                }else{ // chain exists
                    auto& chain_coords = get<0>(chainID_to_coords[chain_ID]);
                    chain_coords.insert(chain_coords.end(), coords.begin(), coords.end());
                }
                progress++;
                cout << progress << "/" << all_intervals.size() << " intervals aligned."<< endl; 
                mtx.unlock();
            } // end of for
        };     // end of the lambda function
        // now distribute the intervals to multiple threads 



        vector<thread> workers;
        for(unsigned int i = 1; i<interval_partition.size(); i++){ // leave one out for the current thread 
            workers.push_back(thread(align_intervals, interval_partition[i]));
        }
        align_intervals(interval_partition[0]);

        for(auto& worker : workers){
            worker.join();
        }

        ofstream ofile_coords;
        ofile_coords.open(coord_ofname);

        for(auto it : chainID_to_coords){
            auto& chainID = it.first;
            auto& coords = get<0>(it.second);
            auto& chain_info = get<1>(it.second); // build1, build1_chr, build2, build2_chr, dir


            // sort by build1 coord
            // within chain, guaranteed to be colinear, so can be sorted using either coordinates 
            sort(coords.begin(), coords.end(),
            [](const tuple<int, int, string, string, string> &a, // build1 coord, build2 coord, Gval, Fval
               const tuple<int, int, string, string, string> &b) -> bool {
                return std::get<0>(a) > std::get<0>(b); 
            });
            int min_coord1 = min(get<0>(coords[0]), get<0>(coords.back())) - window_width/2;
            int max_coord1 = max(get<0>(coords[0]), get<0>(coords.back())) + window_width/2;
            int min_coord2 = min(get<1>(coords[0]), get<1>(coords.back())) - window_width/2;
            int max_coord2 = max(get<1>(coords[0]), get<1>(coords.back())) + window_width/2;
            ofile_coords << ">" << "\t"  << chainID << "\t" << chain_info[0] << "\t" << chain_info[1] << "\t" << min_coord1<<"\t"<< max_coord1 
                                            <<"\t" << chain_info[2] << "\t" << chain_info[3] << "\t" << min_coord2 << "\t" <<max_coord2
                                            <<  "\t" << chain_info[4] << endl;


            for(auto& coord : coords){
                ofile_coords << get<0>(coord) << "\t" << get<1>(coord) << "\t"<<get<2>(coord) << "\t"<< get<3>(coord) << "\t" << get<4>(coord) << endl;
            }
 
            ofile_coords << endl;
        }
        ofile_coords.close();

    } else if (input_file_type == "2") { // positional argument: .bed. extra

        if (cfile_name == "NULL") {
            cout
                << "For input type 1, option -c (build_to_build.coord file) is "
                   "required"
                << endl;
	    help();
            exit(1);
        }



        if (optind + 1 != argc) {
            cout << "One positional argument is required: query bed file"
                 << endl;
	    help();
            exit(1);
        } else {
            ifile_name = argv[optind];
        }

        check_file(ifile_name);
        ifstream ifile(ifile_name);

        if (qbuild == "NULL") {
            cout << "For type -2, must specify the name of the query build (e.g. mm10) after -q" << endl;
	    help();
            exit(1);
        }

        string init_dir = dirname(ifile_name);
        if (ofile_directory ==  "") { // if not specified, save it where the input file is located
            ofile_directory = init_dir;
        }


	if(unique_mapping + multiple_mapping != 1){
	    cout << "Exactly one of -m (multiple) or -u (unique) must be provided" << endl;
	    help();
	    exit(1);
	}

        string ofname_mapped;
        string ofname_nmapped;

	if(ofname_prefix == "NULL"){
		ofname_prefix = basename(ifile_name);
	}
	
        if(multiple_mapping){
            ofname_mapped = ofile_directory + ofname_prefix + ".multiple_mapped";
            ofname_nmapped = ofile_directory + ofname_prefix + ".multiple_not_mapped"; 
        } else if(unique_mapping) {
            ofname_mapped = ofile_directory + ofname_prefix + ".unique_mapped";
            ofname_nmapped = ofile_directory + ofname_prefix + ".unique_not_mapped"; 
        }


        ofstream ofile_mapped;
        ofstream ofile_nmapped;
        ofile_mapped.open(ofname_mapped);
        ofile_nmapped.open(ofname_nmapped);




       // input bed file either has unique IDs for each bed lines in the fourth column, or empty
        // first, read in all the bed lines 
        vector<tuple<string, int, int, string>> bed_inputs;
        vector<string> qchrom_list; // to read in only the chain lines that are necessary
	string line;
        while (!ifile.eof()){
            getline(ifile, line);
            if(ifile.eof()){continue;}

            string chr;
            int start; //  start coordinate for each bed lines
            int end; // end "
	    string bed_id = "";
            istringstream ss(line);
            ss >> chr >> start >> end >> bed_id;

            if(bed_id.length() == 0){ // if no id given, make one 
                bed_id = chr + ":" + to_string(start) + "-" + to_string(end);
            }           
         
            bed_inputs.push_back(make_tuple(chr, start, end, bed_id));
            qchrom_list.push_back(chr);
        }

        cout<<"Initializing Mapper"<<endl;
        Mapper mp(cfile_name, qbuild, qchrom_list);
	cout<<"Mapper loaded" << endl;

        for(auto& query_bed : bed_inputs){
            string chr = get<0>(query_bed);
            int start = get<1>(query_bed);
            int end = get<2>(query_bed); 
            int width = (end-start);

            string bed_id = get<3>(query_bed);
            vector<string> matched_chains = mp.identify_chains(chr, start, end); 

	    // no matching chain 
	    if(matched_chains.size() == 0){
                ofile_nmapped << "# deleted in " << mp.get_tbuild()  <<endl;
                ofile_nmapped << chr + "\t" << start  << "\t" << end << "\t" <<bed_id<<endl;
	        continue;
            }
	    
	    // First, obtain a list of all mappings of the given bed
	    vector<tuple<string, int, int, string>> mapping_list_init = {}; // container for mapped loci
	    // loop through all mapped chains 
            for(unsigned int i = 0; i<matched_chains.size(); i++){
                string chain_id = matched_chains[i];
                vector<tuple<int, int, string>>& coords = mp.get_coords(chain_id);

                tuple<int, int, string> coord = mp.search_query(coords, start, end);
                if(get<0>(coord) != -1){
                    mapping_list_init.push_back(make_tuple(mp.get_tchrom(chain_id), max(0, get<1>(coord) - width/2), get<1>(coord)+ width/2, get<2>(coord)));
                }
            }


	    // scan the initial mapping list, and search for trivial duplicates (significant overlap with other mapping)
	    // sig overlap : overlapped bp is at least half the size of the query bed. 
	    // Merge the trivial duplicates 
	    vector<tuple<string, int, int, string>> mapping_list_final = {}; 
            for(auto& map1 : mapping_list_init){
	        bool trivial_dup = false;
		unsigned int i = 0;
		// check existence for trivial duplication then merge. 
		while(i<mapping_list_final.size() && !trivial_dup){
		    auto& map2 = mapping_list_final[i];
		    bool same_chrom = (get<0>(map1) == get<0>(map2)); 
		    int bp_overlap = min(get<2>(map1), get<2>(map2)) - max(get<1>(map1), get<1>(map2));
		    bool sig_overlap = (bp_overlap > width/2); 
		    if(same_chrom && sig_overlap){
			trivial_dup = true; 
	                get<1>(map2) = min(get<1>(map1), get<1>(map2));
			get<2>(map2) = max(get<2>(map1), get<2>(map2));
		    }
		    i++;
		}


		if(!trivial_dup){
		    mapping_list_final.push_back(map1);
	        }
	    }

            int ncopy = mapping_list_final.size();
	    // matched chain but the input bed line doesn't overlap with any coords
	    if(ncopy == 0){
                ofile_nmapped << "# deleted in " << mp.get_tbuild()  <<endl;
                ofile_nmapped << chr + "\t" << start  << "\t" << end << "\t" <<bed_id<<endl;
                continue;
            }

            if(multiple_mapping){
                for(unsigned int i = 0; i < mapping_list_final.size(); i++){
		    auto& mapped_locus = mapping_list_final[i]; 
                    ofile_mapped << get<0>(mapped_locus) << "\t" << get<1>(mapped_locus)
                                 << "\t" << get<2>(mapped_locus) << "\t" << bed_id
                                 <<  "\t" << i+1 << "\t" << get<3>(mapped_locus)  <<endl;
                }
            } else if(unique_mapping){
	        if(ncopy == 1){
		    auto& mapped_locus = mapping_list_final[0];
		    ofile_mapped << get<0>(mapped_locus) << "\t" << get<1>(mapped_locus)
                                 << "\t" << get<2>(mapped_locus) << "\t" << bed_id
                                 << "\t" << get<3>(mapped_locus)  <<endl;
                }else if (ncopy>=2){ 
                    ofile_nmapped << "# duplicated in " <<  mp.get_tbuild() <<endl;
                    ofile_nmapped << chr << "\t" << start << "\t" << end << "\t" << bed_id<<endl;
                }else{
                    cout << "debug needed (MAPPER-0mapped)" << endl;
		    exit(1);
                }
            } else{
		    cout << "debug required" << endl;
		    exit(1);
	    }
	}
        cout << "Finished mapping" << endl;
        ofile_mapped.close();
        ofile_nmapped.close();
        ifile.close();
    }
    return 0;
}
