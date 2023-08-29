#include "header.h"

// constructor
// current version keeps input bed and out bed width identical. (e.g. 300bp --> 300bp)
//cfile = coordinate chain file 
Mapper::Mapper(string cfile_name, string qbuild, vector<string> qchrom_list) {
    check_file(cfile_name);
    ifstream cfile(cfile_name);

    // format
    // >	1   mm10    chr7    103847958       103886990       hg38    chr11   5267504 5306699 same_strand
    // 103847958       5267504 0.0378858
    // 103847978       5267524 0.0334759

    // chain_id build1 build1_chr range_start range_end build2 build2_chr range_start range_end dir


    // save the alignment result (build_to_build.coord) into a map
    // Partition by chromosome, and sort by position for each chroms.

    // for chain info line
    string trash;
    string chain_id;
    string build1;
    string build1_chr;
    int build1_start;
    int build1_end;
    string build2;
    string build2_chr;
    int build2_start;
    int build2_end;
    string rel_strand;
   

    // for coordinate lines 
    int coord1; 
    int coord2;
    string cons_vals; // for unweighted alignment, one float. For weighted, three floats

    string qchrom; 

    string line;
    bool skip = true; 

    while (getline(cfile, line)) {
        istringstream ss(line);
        if(line[0] == '>'){ // new chain
	    skip = true; 
            ss >> trash >> chain_id >>build1 >> build1_chr >> build1_start >> build1_end
               >> build2 >> build2_chr >> build2_start >> build2_end >> rel_strand;

            if(build1 == qbuild){
                qdim = 0;
                qbuild = build1;
                tbuild = build2;
                qchrom = build1_chr; 
           } else if(build2 == qbuild){
                qdim = 1;
                qbuild = build2;
                tbuild = build1;
                qchrom = build2_chr;
            } else {
                cout << "ERROR: " << qbuild << " does not match any of the builds in the input chain file" << endl;
                exit(1);
            }
 
            if(find(qchrom_list.begin(), qchrom_list.end(), qchrom) != qchrom_list.end()){ // chain lies in query relevant chromosmes
	        skip = false; 
            // fill in chrom_to_cinfo
                if(qdim == 0){
                    if (chrom_to_cinfo.find(build1_chr) == chrom_to_cinfo.end()){
                        chrom_to_cinfo[build1_chr] = {make_tuple(chain_id, build1_start, build1_end)};
                    } else {
                        (chrom_to_cinfo[build1_chr]).push_back(make_tuple(chain_id, build1_start, build1_end));
                    }

                    chainID_to_tchrom[chain_id] = build2_chr;

                } else {
                    if (chrom_to_cinfo.find(build2_chr) == chrom_to_cinfo.end()){
                        chrom_to_cinfo[build2_chr] = {make_tuple(chain_id, build2_start, build2_end)};
                    } else {
                        (chrom_to_cinfo[build2_chr]).push_back(make_tuple(chain_id, build2_start, build2_end));
                    }
                    chainID_to_tchrom[chain_id] = build1_chr;
                }
	    }

        // fill in chainID_to_coords. for coords, first dimension is the query dimension.
        } else {
	    if(!skip){
                ss >> coord1 >> coord2;
	        getline(ss, cons_vals);
                if(qdim == 0){
                    if(chainID_to_coords.find(chain_id) == chainID_to_coords.end()){
                        chainID_to_coords[chain_id] = {make_tuple(coord1, coord2, cons_vals)};
                    } else{
                        (chainID_to_coords[chain_id]).push_back(make_tuple(coord1, coord2, cons_vals)); 
                    }            
                } else {
                    if(chainID_to_coords.find(chain_id) == chainID_to_coords.end()){
                        chainID_to_coords[chain_id] = {make_tuple(coord2, coord1, cons_vals)};
                    } else{
                        (chainID_to_coords[chain_id]).push_back(make_tuple(coord2, coord1, cons_vals));
                    }

                }
	    }
        }
    }

    // now, sort coords for each chain. sort by query coordinate
    for(auto& item : chainID_to_coords){
        auto& coords = item.second;
        sort(coords.begin(), coords.end(),
             [](const tuple<int, int, string> &a, const tuple<int, int, string> &b) -> bool {
                return get<0>(a) < get<0>(b);
        });
        
    }
    cfile.close();
}

// return a list of chains that overlap with the given bed coordinate 
// uses data stored in chrom_to_cinfo
// // chrom -> vector of {chain_id, range_begin, range_end}
vector<string> Mapper::identify_chains(string chr, int start, int end){
    vector<tuple<string, int, int>>& chain_infos = chrom_to_cinfo[chr];

    vector<string> output_chain_list = {};
    for(auto& chain_info : chain_infos){
        bool overlap = min(end, get<2>(chain_info)) > max(start, get<1>(chain_info)); 
        if(overlap){
            output_chain_list.push_back(get<0>(chain_info));
        } 
    }
    return output_chain_list;
}

vector<tuple<int, int, string>>& Mapper::get_coords(string chain_id){
    return chainID_to_coords[chain_id];
}


// input coord is sorted by the first element (query coord)
// obtain all query coordinates  that is contained within the query range 
tuple<int, int, string> Mapper::search_query(vector<tuple<int, int, string>>& coord_list, int qstart, int qend) {
    // start with binary search
    int lower = 0;
    int upper = coord_list.size();
    int x;
    int val;

    int query = (qstart+qend)/2;
    bool found = false;
    while (lower+1 < upper) {
        x = lower + (upper - lower) / 2;
        val = get<0>(coord_list[x]);
        if (query == val) {
            found = true;
            break; // found. But proceed to scan its neighbors to find one with the highest matrix
        } else if (query > val) {
            lower = x;
        } else if (query < val) {
            upper = x;
        }
    }

    int index = -1; // index to the current best matching coordinate to the target genomic coordinate.  
    if(found){ // exact match found
        index = x;
    } else if((qstart < get<0>(coord_list[lower]) && qend > get<0>(coord_list[lower])) ||
            (qstart < get<0>(coord_list[upper]) && qend > get<0>(coord_list[upper])) ){ 
    // lower and upper index differ by one. Either lower or upper has acceptable search solution (contained within query range)
        int d_lower = abs(query - get<0>(coord_list[lower]));
        int d_upper = abs(query - get<0>(coord_list[upper]));
        if(d_lower < d_upper){ // lower is better
            index = lower;
        }else{ // upper is better
            index = upper;
        }
    } else{ // no acceptable match found 
        return make_tuple(-1,-1,"0"); 
    }



    return coord_list[index];

}


string Mapper::get_tchrom(string chainID){
    return chainID_to_tchrom[chainID];
}

string Mapper::get_tbuild(){
    return tbuild;
}
