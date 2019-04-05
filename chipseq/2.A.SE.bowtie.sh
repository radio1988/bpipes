#!/bin/bash
#BSUB -J "bowtie_SE[1-4]"
#BSUB -R rusage[mem=20000]
#BSUB -n 16
#BSUB -q short
#BSUB -W 4:00
##BSUB -N
#BSUB -P "Nathan"

#set -u
#set -e
#set -o pipefail


source config_2018_02_16.sh

cd $webDir

mkdir -p scripts
mkdir -p scripts/log
mkdir -p scripts/log/bowtie2


#BSUB -o /nl/umw_nathan_lawson/pub/2018_02_18_NR2F2_CHIPSEQ/scripts/log/bowtie2/SE.out.%J.%I.txt
#BSUB -e /nl/umw_nathan_lawson/pub/2018_02_18_NR2F2_CHIPSEQ/scripts/log/bowtie2/SE.err.%J.%I.txt


module load bowtie2/2.3.2
module load samtools/1.4.1
module load java/1.8.0_77
module load gnuplot/4.6.5
module load R/3.4.0
module load python3/3.5.0
module load zlib/1.2.8


species=hg19

mkdir -p $ALIGNFOLDER1

# bowtie2 -p 8 -q -t --met-stderr \
#         -x $mm10Bowtie2Index \
#         -X 2000 \
#         -1 fastq/${sample[$i]}.fastq.gz \
#         -2 \
#     |samtools view -bhS -q 20 > bam/${sample[$i]}.bam


# cut and run

bowtie2 --local --very-sensitive-local --no-unal \
-t -q --threads 16 \
-x $hg19Bowtie2Index \
-U $source1/${fq[$i]} \
&2> ${webDir}/${ALIGNFOLDER1}/SE_${sample[$i]}_bowtie2_log.txt \
| samtools view -@ 16 -bhS -q 20 > ${webDir}/${ALIGNFOLDER1}/${sample[$i]}.bam

samtools sort -@ 16 ${webDir}/${ALIGNFOLDER1}/${sample[$i]}.bam -o ${webDir}/${ALIGNFOLDER1}/${sample[$i]}.sort.bam
samtools index ${webDir}/${ALIGNFOLDER1}/${sample[$i]}.sort.bam

if [ $? -eq 0 ]
then
echo "removing the unsorted bam file"
rm ${webDir}/${ALIGNFOLDER1}/${sample[$i]}.bam
else
echo "something wrong with bam sort" >&2
fi

samtools idxstats ${webDir}/${ALIGNFOLDER1}/${sample[$i]}.sort.bam > ${webDir}/${ALIGNFOLDER1}/${sample[$i]}.idxstats.txt
samtools flagstat ${webDir}/${ALIGNFOLDER1}/${sample[$i]}.sort.bam > ${webDir}/${ALIGNFOLDER1}/${sample[$i]}.flagstats.txt



java -Xmx16g -jar $PICARD MarkDuplicates \
I=${webDir}/${ALIGNFOLDER1}/${sample[$i]}.sort.bam \
O= ${webDir}/${ALIGNFOLDER1}/${sample[$i]}.sort.markDup.bam \
M=${webDir}/${ALIGNFOLDER1}/${sample[$i]}.marked_dup_metrics \
REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT

samtools index ${webDir}/${ALIGNFOLDER1}/${sample[$i]}.sort.markDup.bam


#$bamCoverage --bam ${webDir}/${ALIGNFOLDER1}/${sample[$i]}.sort.markDup.bam -o ${webDir}/${ALIGNFOLDER1}/${sample[$i]}.bw \
#--binSize 10 \
#--normalizeUsingRPKM \
#--ignoreForNormalization chrX
#--extendReads

java -Xmx8g -jar $PICARD CollectMultipleMetrics \
I=${webDir}/${ALIGNFOLDER1}/${sample[$i]}.sort.bam  \
O=${webDir}/${ALIGNFOLDER1}/${sample[$i]}.multiple_metrics \
R=$hg19fa
