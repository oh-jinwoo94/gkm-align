#4/21/2020

# K-mean cluster its anchors, where #cluster is determined by the average width and anchor number of the chain.
# k = ceiling(avg_W/L), where L = target avg width of sub-chains for each cluster on average
# After clustering, use medians of the clusters (recall that the anchors are colinear), as "the anchor of the anchors", and generate a list ofr end-to-end coordinates connecting the medians.
# Treat the two anchors at both extremes as medians too.


# genearte 2align for for pairwise end2end alignment between neighboring anchors ()


# .2align file for mat
#species 1 = query. species 2 = target 
#species_1_genome_build  species_1_chr   species_1_range_begin   species_1_range_end     species_2_genome_build  species_2_chr   species_2_range_begin   species_2_range_end     relative_strand description chain_ID
#hg38    chr9    98248109        98328109        mm10    chr4    46622305        46702305        same_strand     ENSG00000136928_ENSMUSG00000039809 chain10

import sys
from sklearn.cluster import KMeans
import numpy as np
import math

def main(argv = sys.argv):
    if(len(argv) != 6):
        print("Usage: {0} {chain file} {L} {E1: chain corner extend} {E2: center corner extend} {ofile name}")
        sys.exit()

    chain_file = open(argv[1], 'r')

    L = int(argv[2]) # avg chain width per centroid
    E1 = eval(argv[3]) 
    E2 = eval(argv[4]) 

    ofile = open(argv[5], 'w')
    chain_dict = {}  # dictionary of  lists of genomic coordinates for anchors. key = chain_info
    chain_info_list = []

    chain_file.readline() # first line is a header
    for i,line in enumerate(chain_file):
        words = line.split()
        tcoord, qcoord, chain_info = int(words[0]), int(words[1]), words[2] 
        if(chain_info not in chain_dict.keys()):
            chain_dict[chain_info] = [[tcoord, qcoord]]
            chain_info_list.append(chain_info)
        else:
            (chain_dict[chain_info]).append([tcoord, qcoord])
            
    ofile.write("\t".join(["target_genome_build", "target_chr", "target_range_begin", "target_range_end", \
                          "query_genome_build",  "query_chr", "query_range_begin", "query_range_end", \
                          "relative_strand", "chainID"]) + '\n')
    # loop through chain. For each chain, get "anchor of the anchors" for end to end alignment (including the leftmost and rightmost extremes) 
    for chain_info in chain_info_list:
        chain = chain_dict[chain_info]
        #target_coordinate       query_coordinate       id 
        #4407547        93723   1

        ID_comp = chain_info.split("/")
        chain_ID = ID_comp[-1]
        tinfo, qinfo = (ID_comp[0]).split("_"), (ID_comp[1]).split("_")

        tbuild, qbuild = tinfo[0], qinfo[0]
        tchr, qchr = "_".join(tinfo[1:]), "_".join(qinfo[1:])

        qStrand = ID_comp[-2]
        tmin, tmax = min(chain, key = lambda x: x[0])[0], max(chain, key = lambda x: x[0])[0]
        qmin, qmax = min(chain, key = lambda x: x[1])[1], max(chain, key = lambda x: x[1])[1]

        t_width, q_width = tmax - tmin, qmax - qmin
        avg_chain_width = (t_width + q_width) / 2
        num_anchor = len(chain)

        num_cluster = int(min(max(math.ceil(avg_chain_width/L-1), 1), num_anchor)) # keep 1<=ncluster<=num_anchor
        X = np.array(chain) # array of 2D coordinates
        kmeans = KMeans(n_clusters=num_cluster, random_state=0).fit(X)
        lab = kmeans.labels_
    

        # generate "anchor of the anchors" using kmeans clustering 
        center_list = []
        for k in range(0, num_cluster):
            # anchors are colinear anyways, and also in the same chromosome Use target genome coordinate to find median.
            X_k = X[lab == k]
            mid = (np.sort(X_k[:,0])[len(X_k[:,0])//2]) # middle element within a cluster
            mid_index = ((list(X[:,0])).index(mid))  # index in the list of all anchor list in a chain
            center_list.append(chain[mid_index]) # this can be done because we know that anchors have unique coordinates within a chain

        # add the two extremes as well, only if num_anchor>1
        # extend by E
        if(num_anchor > 1):
            extreme_1, extreme_2 = list(chain[0]), list(chain[-1])
            if(qStrand == "+"):
                extreme_1[0] -= E1
                extreme_1[1] -= E1
            else:
                extreme_1[0] += E1
                extreme_1[1] -= E1
            center_list.append(extreme_1)

            if(qStrand == "+"):
                extreme_2[0] += E1
                extreme_2[1] += E1
            else:
                extreme_2[0] -= E1
                extreme_2[1] += E1
            center_list.append(extreme_2)

        center_list = sorted(center_list, key=lambda center: center[0]) # sort by target coordinate
        for i in range(0, len(center_list)-1):
            center_1, center_2 = center_list[i], center_list[i+1]
            tcoord_1, qcoord_1 = center_1[0], center_1[1]
            tcoord_2, qcoord_2 = center_2[0], center_2[1]
            if(qStrand == "+"):
                ofile.write("\t".join([tbuild, tchr, str(tcoord_1-E2), str(tcoord_2+E2), qbuild, qchr, str(qcoord_1-E2), str(qcoord_2+E2), \
                                 "same_strand",  chain_ID]) + '\n') 
            else:
                # for "-" query coordinate decreases with anchor order 
                ofile.write("\t".join([tbuild, tchr, str(tcoord_1-E2), str(tcoord_2+E2), qbuild, qchr, str(qcoord_2-E2), str(qcoord_1+E2), \
                                 "diff_strand",  chain_ID]) + '\n')
        
        # loop through all the intervals (with overlapping E extension)
        

    ofile.close()
    chain_file.close()
main()
