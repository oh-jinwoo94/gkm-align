seq1 = ""
with open("test_hg38_gkm-masked_h.fa", 'r') as ifile:
    for i, line in enumerate(ifile):
        if(i%2 == 1):
            seq1 += line.rstrip()


seq2 = ""
with open("test_hg38_gkm-masked_m.fa", 'r') as ifile:
    for i, line in enumerate(ifile):
        if(i%2 == 1):
            seq2 += line.rstrip()


right = 0
for i in range(0, len(seq1)):
    if(seq1[i] == seq2[i]):
        right += 1

print(f'{round(right/len(seq1)*100, 2)}% match')

