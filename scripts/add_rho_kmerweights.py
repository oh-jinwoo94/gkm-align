# Given weight file and model file, add rho in model file to all kmer weights.
#nr11mers_genomic_model_mm10_l300_n30000_r10_weights.out
#genomic_model_mm10_l300_n30000_r10.model.txt


import sys

def main(argv = sys.argv):
    if(len(argv) != 4):
        print("<0> <weight file> <model.txt> <ofile name>")
        sys.exit(1)

    wfile = open(argv[1], 'r')
    mfile = open(argv[2], 'r')
    ofile = open(argv[3], 'w')

    rho = "N/A"
    for line in mfile:
        words = line.split()
        if(words[0] == "rho"):
            rho = eval(words[1])
            break
    if(rho == "N/A"):
        print("incorrect file format")
        sys.exit(1)
     

    for line in wfile:
        words = line.split()
        ofile.write("\t".join([words[0], str(eval(words[1]) + rho)]) + '\n')

    ofile.close()
    mfile.close()
    wfile.close()
main()
