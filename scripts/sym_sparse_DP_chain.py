#!/usr/bin/env python3.6

# 4/28/2020
# symmetric version 
#Have two sorted lists of anchors sorted by both query or target location.  At each iteration, consider the leftmost anchor in the sorted lists. If the two anchors are predicted to have the best score in the same chain, use the one that has higher eventual score. Then, remove that anchor from both of the lists. Continue. 



import math
import sys

class Anchor:
    def __init__(self, t, q): # t: target genome loc q: query genome loc
        self.t=t
        self.q=q
    def __eq__(self, other):
        if(self.t == other.t and self.q == other.q):
            return True
        else:
            return False  

    def toString(self):
        return str(self.t) + "_" + str(self.q)

class Anchor_Set: 
#a group list of anchors that have the same orientation & belong to same target AND query chromosomes.
    def __init__(self, anchor, orientation):
        self.anchor_list = [anchor]
        self.orientation = orientation
    def add_anchor(self, anchor):
        (self.anchor_list).append(anchor)

# a list of anchors. Its location defined as the location of its head (the latest addition to the chain)
class Chain:
    def __init__(self, anchor):
        self.anchor_list = [anchor]
        self.t, self.q = anchor.t, anchor.q
        self.score = 1

    def add_anchor(self, anchor, score): # add a new anchor, and set it as the new head
        if(anchor.q<self.q): # not colinear
            print("ERROR: DEBUGGING REQUIRED.")
            sys.exit()
        self.anchor_list.append(anchor)
        self.t, self.q = anchor.t, anchor.q
        self.score = score 

    def find_width_ratio(self):

        if(len(self.anchor_list) == 1):
            return 1
        else:
            first, last = self.anchor_list[0], self.anchor_list[-1]
            w1 = abs(first.q-last.q)
            w2 = abs(first.t-last.t)
            if(w1>w2):
                return float(w2) / float(w1)
            else:
                return float(w1) / float(w2)


def distance(A, B):
   return(abs(A.t-B.t) + abs(A.q-B.q)) 


# returns index for the best chain and the resultant score. -1, -1 if none (self chain)
def find_best_chain(anchor, orientation, chain_list, indel, max_dist):
    max_score = 1 # score of making its own chain
    best_chain = -1
    for j, chain in enumerate(chain_list):
        colinear = (orientation=="+" and (anchor.t > chain.t and anchor.q > chain.q)) or \
                   (orientation=="-" and (anchor.t < chain.t and anchor.q > chain.q))

        close_enough = (distance(chain, anchor) < max_dist)
        if(colinear and close_enough):
            pass #viable chain
        else:
            continue #skip this chain

        score = chain.score + 1 - indel*distance(chain, anchor)
        if(score > max_score):
            max_score = score
            best_chain = j
#        if(orientation=="-" and max_score>0):
#            print(str(chain.q)+"_"+str(anchor.q)) 
    return best_chain, max_score


def find_chains(anchor_set, indel, max_dist):
    chain_list = []
    completed = {anchor.toString() : 0 for anchor in anchor_set.anchor_list} # contains information on whether an anchor has been added to a chain

    qsorted_anchor_list = sorted(anchor_set.anchor_list, key = lambda x:x.q)
    if(anchor_set.orientation == "+"):
        tsorted_anchor_list = sorted(anchor_set.anchor_list, key = lambda x:x.t)
    else: 
        tsorted_anchor_list = sorted(anchor_set.anchor_list, key = lambda x:(-x.t))
    q_index, t_index = 0, 0 

    # start with the first anchor in query (this level of asymmetry almost has no effect)
    qanchor = qsorted_anchor_list[q_index]
    chain_list.append(Chain(qanchor))
    completed[qanchor.toString()] = 1
    q_index += 1


    N = len(anchor_set.anchor_list)
    # iterate until one of the sorted lists reaches its end
    while(q_index<N and t_index<N):

        qanchor, tanchor = qsorted_anchor_list[q_index], tsorted_anchor_list[t_index] 

        if(qanchor == tanchor): # if the two indices point to the same anchor, do it once for qanchor 
            if(completed[qanchor.toString()] == 0): # proceed if this anchor has not been encountered yet
                q_best_chain, q_max_score = find_best_chain(qanchor, anchor_set.orientation, chain_list, indel, max_dist)
                if(q_best_chain == -1):
                    chain_list.append(Chain(qanchor))
                    completed[qanchor.toString()] = 1 
                else:
                    (chain_list[q_best_chain]).add_anchor(qanchor, q_max_score)
                    completed[qanchor.toString()] = 1
            q_index+=1
            t_index+=1 

        else:   # if point to different anchors
        
            # skip anchor if already added (i.e. from the other sorted list)
            if(completed[qanchor.toString()] == 1):
                q_index+=1

            elif(completed[tanchor.toString()] == 1):
                t_index+=1
 
            else: #if the anchors pointed by qanchor and tanchor both have not been encountered
                q_best_chain, q_max_score = find_best_chain(qanchor, anchor_set.orientation, chain_list, indel, max_dist)
                t_best_chain, t_max_score = find_best_chain(tanchor, anchor_set.orientation, chain_list, indel, max_dist)
    
                # if the best chains for both of the anchors are same,
                if(q_best_chain == t_best_chain): 
                    if(q_best_chain == -1): # if they both have higher scores as its solo chain
                        chain_list.append(Chain(qanchor))
                        completed[qanchor.toString()] = 1
                        chain_list.append(Chain(tanchor))
                        completed[tanchor.toString()] = 1
                        q_index += 1
                        t_index += 1

                    elif(q_max_score >= t_max_score): # if adding query anchor has higher score
                        (chain_list[q_best_chain]).add_anchor(qanchor, q_max_score)
                        completed[qanchor.toString()] = 1
                        q_index += 1
                    else: # if adding tanchor has higher score
                        (chain_list[t_best_chain]).add_anchor(tanchor, t_max_score)
                        completed[tanchor.toString()] = 1
                        t_index += 1
                else: # if the two anchors are headed to different chains, add them independently 

                    if(q_best_chain == -1): # if its has higher score as a solo chain.
                        chain_list.append(Chain(qanchor))
                        completed[qanchor.toString()] = 1
                    else:
                        (chain_list[q_best_chain]).add_anchor(qanchor, q_max_score)
                        completed[qanchor.toString()] = 1

                    if(t_best_chain == -1):
                        chain_list.append(Chain(tanchor))
                        completed[tanchor.toString()] = 1
                    else:
                        (chain_list[t_best_chain]).add_anchor(tanchor, t_max_score)
                        completed[tanchor.toString()] = 1

                    q_index += 1 
                    t_index += 1
    return chain_list



def main(argv = sys.argv):

    if(len(argv) != 11):
        print("Usage: {0} {anchor file} {target genome id} {query genome id} {target genome chrom.sizes} {query genome chrom.size} {indel penalty} {maximum distance} {ratio threshold} {size threshold}  {ofile name}")
        sys.exit()

    # anchor file 
    ifile = open(argv[1], 'r')
    
    tgenome, qgenome = argv[2], argv[3]

    # read target chrom sizes
    t_chrom_size = {}
    with open(argv[4]) as sz_file:
        for line in sz_file:
            words = line.split()
            t_chrom_size[words[0]] = int(words[1]) 

    # read query chrom sizes
    q_chrom_size = {}
    with open(argv[5]) as sz_file:
        for line in sz_file:
            words = line.split()
            q_chrom_size[words[0]] = int(words[1])

    indel_penalty = eval(argv[6])
    max_dist = eval(argv[7])
    rthresh = eval(argv[8])
    sthresh = eval(argv[9])

    ofile = open(argv[10], 'w')


    # keys: chromosome pairs with strand info. e.g. chr8_chr10_+
    anchor_partition = dict()
    
    for line in ifile:
        #9 chr1 174089066 174089450 chr1 62973 63357 + 5496
        words = line.split()
        if(line[0] == "#" or len(words) <= 1): # only need seq locs
            continue
        chrom_t, chrom_q, orientation = words[1], words[4], words[7]

        if((chrom_t not in t_chrom_size.keys()) or (chrom_q not in q_chrom_size.keys())): # skip over cryptic chromosome ids 
            continue
        if(orientation == '+'):
            tloc = (int(words[2]) + int(words[3])) // 2
            qloc = (int(words[5]) + int(words[6])) // 2
        elif(orientation =='-'):
            tloc = (int(words[2]) + int(words[3])) // 2
            qloc = q_chrom_size[chrom_q] - (int(words[5]) + int(words[6])) // 2 - 1
        #print(str(tloc) + "\t" + str(qloc))
        part_id = "/".join([tgenome + "_" + words[1], qgenome + "_" + words[4], orientation])
        if(part_id not in anchor_partition.keys()):
            anchor_partition[part_id] = Anchor_Set(Anchor(tloc, qloc), orientation)
        else:
            (anchor_partition[part_id]).add_anchor(Anchor(tloc, qloc))
  
    ofile.write("\t".join(["target_coordinate", "query_coordinate", "chain_ID"]) + '\n')
    chain_id = 1 
    for part_id in anchor_partition.keys():
        anchor_set = anchor_partition[part_id]
        chains = find_chains(anchor_set, indel_penalty, max_dist)  
        for chain in chains:
            # large enough and width height ratio is reasonable 
            if((len(chain.anchor_list) >= sthresh) and (chain.find_width_ratio() >= rthresh)):
                for node in chain.anchor_list:
                    chain_info = "/".join([part_id, "chain" + str(chain_id)])
                    ofile.write("\t".join([str(node.t), str(node.q), str(chain_info)]) + '\n')
            chain_id = chain_id+1
    ifile.close()
    ofile.close()
if __name__=="__main__":
    main()
