#include "header.h"

// constructor
// current version keeps input bed and out bed width identical. (e.g. 300bp --> 300bp)
//cfile = coordinate chain file 
Mapper::Mapper(string cfile_name, string qbuild, vector<string> qchrom_list) {
    check_file(cfile_name);
    ifstream cfile(cfile_name);

    string line;
    
    // cache - pointer to the vector we are currently filling.
    vector<tuple<int, int, string>>* current_coords = nullptr;

    // temporary variables for parsing header
    string trash, chain_id, build1, build1_chr, build2, build2_chr, rel_strand;
    int build1_start, build1_end, build2_start, build2_end;
    string qchrom;
    
    // Parsing state
    bool skip = true;

    while (getline(cfile, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') { // header line 
            istringstream ss(line);
            ss >> trash >> chain_id >> build1 >> build1_chr >> build1_start >> build1_end
               >> build2 >> build2_chr >> build2_start >> build2_end >> rel_strand;

            if (build1 == qbuild) {
                qdim = 0; qbuild = build1; tbuild = build2; qchrom = build1_chr;
            } else if (build2 == qbuild) {
                qdim = 1; qbuild = build2; tbuild = build1; qchrom = build2_chr;
            } else {
                cout << "ERROR: " << qbuild << " does not match any of the builds in the input chain file" << endl;
                exit(1);
            }

            if (find(qchrom_list.begin(), qchrom_list.end(), qchrom) != qchrom_list.end()) {
                skip = false;
                
                // Update Metadata
                if (qdim == 0) {
                    chrom_to_cinfo[build1_chr].emplace_back(chain_id, build1_start, build1_end);
                    chainID_to_tchrom[chain_id] = build2_chr;
                } else {
                    chrom_to_cinfo[build2_chr].emplace_back(chain_id, build2_start, build2_end);
                    chainID_to_tchrom[chain_id] = build1_chr;
                }

                // POINTER CACHE: Set the pointer to the current vector
                current_coords = &chainID_to_coords[chain_id];
                
            } else {
                skip = true;
                current_coords = nullptr; 
            }

        } else { // coordinate lines
            if (!skip && current_coords != nullptr) {
                const char* p = line.c_str();
                char* end_ptr;

                // fast integer parsing
                int c1 = strtol(p, &end_ptr, 10);
                p = end_ptr;
                int c2 = strtol(p, &end_ptr, 10);
                
                // fast string capture (rest of line)
                while (*end_ptr && isspace(*end_ptr)) end_ptr++;
                string cons_vals(end_ptr); 

                if (qdim == 0) {
                    current_coords->emplace_back(c1, c2, cons_vals);
                } else {
                    current_coords->emplace_back(c2, c1, cons_vals);
                }
            }
        }
    }
    
    // sort logic
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
