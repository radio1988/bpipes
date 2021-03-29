from Bio import SeqIO
import sys

if ( len(sys.argv) != 4 ):
    print ('usage: fasta_padding.py <input.fa> <target-length> <output.fa>')
    print ('e.g.: fasta_padding.py input.fa 100 output.fa')
    print('cmd:', sys.argv)
    exit("wrong cmd")
else:
    print('cmd:', sys.argv)

LT = int(sys.argv[2]) # length target

OUT = open(sys.argv[3],'w')

for record in SeqIO.parse(sys.argv[1], "fasta"):
    L = len(record.seq)
    seq = record.seq
    out = ""
    if L >= LT:
        seq = record.seq[(L//2 - LT//2):(L//2 + LT//2)]
        out = seq.upper()
    else:
        head = "N" * ((LT - L) // 2)
        out = head + seq
        tail = "N" * (LT - len(out))
        out = head + seq + tail
        out = out.upper()
    outstr = ">" + record.id + "\n" + out + "\n"
    outstr = str(outstr)
    OUT.write(outstr)
