import pysam
from Bio import SeqIO
import sys
import re
import os

# need source activate damid
print("usage: \
    python filter_bam.py <sample.bam> <genome.fa> <motif> <out-bam-name>")
print("example: \
    python filter_bam.py sorted_reads/L7_S2.bam hg38.primary.fa \
    GATC DamID_filtered_reads/L7_S2.bam")
print("command:", sys.argv)
outname=sys.argv[4]
print("input:", sys.argv[1])
print("output:", outname)
print("caveat: motif must be palindrome")
print("filter_bam.py will remove reads with \
    MQ < 20, secondary alignments, \
    and reads without motif after cut-site on the genome")
print("input bam must be sorted")

def contain_motif(motif="GATC", seq="AcTggatc"):
    motif = motif.upper()
    seq = seq.upper()
    if motif in seq:
        return True
    else:
        return False    

# read genome fasta
genome_dict = {rec.id : rec.seq for rec in SeqIO.parse(sys.argv[2], "fasta")}
#chrom = str(seq_dict['chr1'])


# go through rawbam, filter, then output fltbam
rawbam = pysam.AlignmentFile(sys.argv[1], "rb")
fltbam =  pysam.AlignmentFile(outname, "wb", template=rawbam)


iter = rawbam.fetch()
n1 = 0
n2 = 0
for x in iter:
    if x.is_secondary:
        continue
    if x.is_unmapped:
        continue
    if x.mapping_quality < 20:
        continue
    n1 += 1
    
    if x.reference_name in genome_dict:
        chrom = genome_dict[x.reference_name]
        if x.is_reverse:
            pos = x.reference_end
            seq = chrom[pos:pos+11] # RC same as Seq for DpnII
        else:
            pos = x.reference_start
            seq = chrom[pos-11:pos]
        
        if contain_motif(sys.argv[3], seq):
            n2 += 1
            fltbam.write(x)  
    else:
        print(x.reference_name, file=sys.stderr)

print("n1 = ", n1, ", n2 = ", n2, "n2/n1", n2/n1)
rawbam.close()
fltbam.close()

print("creating index...")
cmd="samtools index " + outname
os.system(cmd)
