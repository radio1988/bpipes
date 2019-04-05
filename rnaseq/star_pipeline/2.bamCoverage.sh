#!/bin/bash
#BSUB -P 'LiR'
#BSUB -J 'bw[1-6]'
#BSUB -R rusage[mem=2000]
#BSUB -n 4
#BSUB -R "span[hosts=1]"
mkdir -p log
#BSUB -o ./log/bamCoverage.%J.%I.log
#BSUB -W 4:00
#BSUB -q short
#BSUB -N
$bam_dir=star

bam_files=(`ls ../$bam_dir/*sort*bam`)
i=$(($LSB_JOBINDEX- 1))

file=${bam_files[$i]}
echo "For: ", $file

module load samtools/dev-2016_06_06
samtools index $file ${file}.bai

module purge
bamCoverage --bam $file -o ${file}.bw \
--numberOfProcessors 4 \
--normalizeUsing RPKM --outFileFormat bigwig 
