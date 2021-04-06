import sys
import math
from Bio import SeqIO



if len(sys.argv) != 5:
    print("usage: python scripts/get_summit_neighbour.py <genome> <input.summit.bed> <width-around_summit> <outout.fa>")
    print("e.g.: python scripts/get_summit_neighbour.py /project/umw_mccb/genome/Homo_sapiens/ucsc_hg38_primary/hg38.primary.fa macs2_DamID_contrast/contrast1/G1_vs_ctrl_summits.bed 1000 macs2_DamID_contrast/contrast1/G1_vs_ctrl_summits.1000.fa")
    exit("cmd error")
else:
    print("cmd:", sys.argv)




genome = sys.argv[1]
summit = sys.argv[2]
width = int(sys.argv[3])
output = sys.argv[4]

# read genome
id2seq = {}
for record in SeqIO.parse(genome, "fasta"):
    id2seq[record.id] = record.seq


# summit file
OUT = open(output, 'w') 
with open(summit , 'r') as f:
    for line in f:
        ele = line.split("\t")
        chr = ele[0]
        location = int(ele[1])
        peak_name = ele[3]
        chr_len = len(id2seq[chr])
        left = max(0, location - math.floor(width/2))
        right = min(location + math.ceil(width/2), chr_len)
        extract = id2seq[chr][left:right]
        if location - width//2 < 0:
            extract = 'N' * (width//2 - location) + extract 
        if location + math.ceil(width/2) > chr_len:
            extract = extract + 'N' * (location + math.ceil(width/2) - chr_len)
        extract = extract.upper()
        extract ='>' + peak_name + "\n" + extract + "\n"
        extract = str(extract)
        OUT.write(extract)
OUT.close()

