import sys
import math
from Bio import SeqIO



if len(sys.argv) != 4:
    print("usage: python scripts/get_peak_fasta.py <genome> <input_clean.narrowPeak> <outout.fa>")
    print("e.g.: python scripts/get_peak_fasta.py genome.fa contrast.narrowPeak  peak.fa")
    exit("cmd error")
else:
    print("cmd:", sys.argv)

genome = sys.argv[1]
peak = sys.argv[2]
output = sys.argv[3]

# read genome
id2seq = {}
for record in SeqIO.parse(genome, "fasta"):
    id2seq[record.id] = record.seq


# peak file
OUT = open(output, 'w') 
with open(peak , 'r') as f:
    for line in f:
        ele = line.split("\t")
        chr = ele[0]
        start = int(ele[1])
        end = int(ele[2])
        peak_name = ele[3]
        chr_len = len(id2seq[chr])
        left = start
        right = end
        extract = id2seq[chr][left:right]
        if start < 0:
            extract = 'N' * (0 - start) + extract 
        if end > chr_len:
            extract = extract + 'N' * (end - chr_len)
        extract = extract.upper()
        extract ='>' + chr + ":" + str(left) + '-' + str(right) + '; ' + peak_name + "\n" + extract + "\n"
        extract = str(extract)
        OUT.write(extract)
OUT.close()

