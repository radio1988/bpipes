#They are located in the following directories (file names listed underneath directory name):
#
#/project/umw_nathan_lawson/deep_seq_data/PreProcess/06FEB18_MiSeq_BLRLC
#TRC-NS-INPUT-LIB1AD1_S1_L001_R1_001.fastq.gz
#TRC-NS-NR2F2-LIB1AD2_S2_L001_R1_001.fastq.gz
#TRC-NS-NR2F2-LIB2AD3_S3_L001_R1_001.fastq.gz
#TRC-NS-NR2F2-LIB3AD4_S4_L001_R1_001.fastq.gz
#
#
#/project/umw_nathan_lawson/deep_seq_data/PreProcess/09FEB18_HS4K
#shVESP-1_input-LIB1_S6_L005_R1_001.fastq.gz
#shVESP-1_input-LIB1_S6_L005_R2_001.fastq.gz
#shVESP-1_NR2F2_IP-LIB1_S7_L005_R1_001.fastq.gz
#shVESP-1_NR2F2_IP-LIB1_S7_L005_R2_001.fastq.gz
#shVESP-1_NR2F2_IP-LIB2_S8_L005_R1_001.fastq.gz
#shVESP-1_NR2F2_IP-LIB2_S8_L005_R2_001.fastq.gz
#shVESP-1_NR2F2_IP-LIB3_S9_L005_R1_001.fastq.gz
#shVESP-1_NR2F2_IP-LIB3_S9_L005_R2_001.fastq.gz
#
#The first set is single end read, the second set is paired end.
#
#Please map onto hg19, then process with MACS (whatever version you are currently using) to identify binding sites
#(input samples are listed - there is only one versus three different ChIPs, but just use this for all for now).
#Also, please provide bigwig files of mapped reads that we can upload to UCSC.
#
#You can place all output files in:  /project/umw_nathan_lawson/lawsonn/pub/NR2F2vespKDchipSeq
#


workingDir='/project/umw_nathan_lawson/Jun/2018_02_16_NR2F2_CHIPSEQ'
outDir="/project/umw_nathan_lawson/lawsonn/pub/NR2F2vespKDchipSeq"
webDir="/nl/umw_nathan_lawson/pub/2018_02_18_NR2F2_CHIPSEQ"

sample=(INPUT NR2F2_1 NR2F2_2 NR2F2_3)

PICARD=/project/umw_mccb/software/picard2.10.0/picard.jar
bamCoverage=/home/jy91w/.local/bin/bamCoverage
hg19fa=/share/data/umw_biocore/genome_data/human/hg19/hg19.fa
hg19Bowtie2Index=/share/data/umw_biocore/genome_data/human/hg19/hg19


source1=/project/umw_nathan_lawson/deep_seq_data/PreProcess/06FEB18_MiSeq_BLRLC

fq=(TRC-NS-INPUT-LIB1AD1_S1_L001_R1_001.fastq.gz \
TRC-NS-NR2F2-LIB1AD2_S2_L001_R1_001.fastq.gz \
TRC-NS-NR2F2-LIB2AD3_S3_L001_R1_001.fastq.gz \
TRC-NS-NR2F2-LIB3AD4_S4_L001_R1_001.fastq.gz)

source2="/project/umw_nathan_lawson/deep_seq_data/PreProcess/09FEB18_HS4K"
fq1=(shVESP-1_input-LIB1_S6_L005_R1_001.fastq.gz \
shVESP-1_NR2F2_IP-LIB1_S7_L005_R1_001.fastq.gz \
shVESP-1_NR2F2_IP-LIB2_S8_L005_R1_001.fastq.gz \
shVESP-1_NR2F2_IP-LIB3_S9_L005_R1_001.fastq.gz)

fq2=(shVESP-1_input-LIB1_S6_L005_R2_001.fastq.gz \
shVESP-1_NR2F2_IP-LIB1_S7_L005_R2_001.fastq.gz \
shVESP-1_NR2F2_IP-LIB2_S8_L005_R2_001.fastq.gz \
shVESP-1_NR2F2_IP-LIB3_S9_L005_R2_001.fastq.gz)

ALIGNFOLDER2=bowtie2_PE
ALIGNFOLDER1=bowtie2_SE

i=$(($LSB_JOBINDEX- 1))
