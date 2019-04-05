#!/bin/bash
#BSUB -J "bowtie[1-3]"
#BSUB -R rusage[mem=16000]
#BSUB -n 8
#BSUB -q long
#BSUB -W 24:00
##BSUB -N
#BSUB -P "Nathan"
#BSUB -R span[hosts=1]
#set -u
#set -e
#set -o pipefail


source config_2018_04_25.sh

cd $workingDir

mkdir -p scripts
mkdir -p scripts/log
mkdir -p scripts/log/bowtie2


#BSUB -o /project/umw_cole_haynes/Jun/Qi/2018_04_25_CHIPSEQ/scripts/log/bowtie2/out.%J.%I.txt
#BSUB -e /project/umw_cole_haynes/Jun/Qi/2018_04_25_CHIPSEQ/scripts/log/bowtie2/err.%J.%I.txt


module load bowtie2/2.3.2
module load samtools/1.4.1
module load java/1.8.0_77
module load gnuplot/4.6.5
module load R/3.4.0
module load python3/3.5.0
module load zlib/1.2.8


species=cel

mkdir -p $ALIGNFOLDER2

# bowtie2 -p 8 -q -t --met-stderr \
# 		-x $mm10Bowtie2Index \
#         -X 2000 \
# 		-1 fastq/${sample[$i]}.fastq.gz \
# 		-2 \
#     |samtools view -bhS -q 20 > bam/${sample[$i]}.bam


# cut and run

bowtie2  --very-sensitive-local --no-unal --no-mixed --no-discordant \
-t -q -I 10 -X 700 --threads 8 \
-x $Bowtie2Index \
-1 $source2/${fq1[$i]} \
-2 $source2/${fq2[$i]} \
2> ${webDir}/${ALIGNFOLDER2}/${sample[$i]}.bowtie2.log.txt \
| samtools view -@ 8 -bhS -q 20 > ${webDir}/${ALIGNFOLDER2}/${sample[$i]}.bam

samtools sort -@ 8 ${webDir}/${ALIGNFOLDER2}/${sample[$i]}.bam -o ${webDir}/${ALIGNFOLDER2}/${sample[$i]}.sort.bam
samtools index ${webDir}/${ALIGNFOLDER2}/${sample[$i]}.sort.bam

if [ $? -eq 0 ]
then
echo "removing the unsorted bam file"
rm ${webDir}/${ALIGNFOLDER2}/${sample[$i]}.bam
else
echo "something wrong with bam sort" >&2
fi

samtools idxstats ${webDir}/${ALIGNFOLDER2}/${sample[$i]}.sort.bam > ${webDir}/${ALIGNFOLDER2}/${sample[$i]}.idxstats.txt
samtools flagstat ${webDir}/${ALIGNFOLDER2}/${sample[$i]}.sort.bam > ${webDir}/${ALIGNFOLDER2}/${sample[$i]}.flagstats.txt



java -Xmx16g -jar $PICARD MarkDuplicates \
I=${webDir}/${ALIGNFOLDER2}/${sample[$i]}.sort.bam \
O= ${webDir}/${ALIGNFOLDER2}/${sample[$i]}.sort.markDup.bam \
M=${webDir}/${ALIGNFOLDER2}/${sample[$i]}.marked_dup_metrics \
REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT

samtools index ${webDir}/${ALIGNFOLDER2}/${sample[$i]}.sort.markDup.bam

mkdir -p ${webDir}
$bamCoverage --bam ${webDir}/${ALIGNFOLDER2}/${sample[$i]}.sort.markDup.bam -o ${webDir}/${ALIGNFOLDER2}/${sample[$i]}.bw \
--binSize 10 \
--normalizeUsingRPKM \
--ignoreForNormalization chrX \
--extendReads

java -Xmx16g -jar $PICARD CollectInsertSizeMetrics \
I=${webDir}/${ALIGNFOLDER2}/${sample[$i]}.sort.bam  \
O=${webDir}/${ALIGNFOLDER2}/${sample[$i]}_insert.size_matrix.txt \
H=${webDir}/${ALIGNFOLDER2}/${sample[$i]}_insert_size_histogram.pdf \
M=0.5

java -Xmx16g -jar $PICARD CollectMultipleMetrics \
I=${webDir}/${ALIGNFOLDER2}/${sample[$i]}.sort.bam  \
O=${webDir}/${ALIGNFOLDER2}/${sample[$i]}.multiple_metrics \
R=$fa

