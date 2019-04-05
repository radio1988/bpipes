#!/bin/bash
#BSUB -J "4.PFinger"
#BSUB -R rusage[mem=20000]
#BSUB -n 16
#BSUB -q short
#BSUB -W 4:00
##BSUB -N
#BSUB -P "jun"
#BSUB -R "span[hosts=1]" # All hosts on the same chassis
####BSUB -w 'ended("bowtie")'


### variables here
bamCoverage=/home/jy91w/.local/bin/bamCoverage
plotFingerprint=/home/jy91w/.local/bin/plotFingerprint
ALIGNFOLDER=bam
OUT=deeptools


workingDir='/project/umw_leslie_shaw/Jun_Allele/Mag/20170720-CUTandRUN'

cd $workingDir

mkdir -p scripts
mkdir -p scripts/log
mkdir -p scripts/log/dpt


#BSUB -o /project/umw_leslie_shaw/Jun_Allele/Mag/20170720-CUTandRUN/scripts/log/pdt/4a.out.%J.%I.txt
#BSUB -e /project/umw_leslie_shaw/Jun_Allele/Mag/20170720-CUTandRUN/scripts/log/dpt/4a.err.%J.%I.txt


module load python3/3.5.0
#module load python3/3.5.0_packages/cython/0.23.5
#module load python3/3.5.0_packages/numpy/1.10.1

#i=$(($LSB_JOBINDEX - 1))

#### sample=(Zfish-10 Zfish-15)

#$bamCoverage --bam bam/${sample[$i]}.sub.sort.markDup.bam -o $dest/${sample[$i]}.bw \
#--binSize 10 \
#--normalizeUsingRPKM \
#--ignoreForNormalization chrX \
#--extendReads

sample=(pre_input pre_IP post_input post_IP)


mkdir -p $OUT

$plotFingerprint \
-b ${ALIGNFOLDER}/pre_input.sort.markDup.bam ${ALIGNFOLDER}/pre_IP.sort.markDup.bam ${ALIGNFOLDER}/post_input.sort.markDup.bam ${ALIGNFOLDER}/post_IP.sort.markDup.bam  \
--labels pre_input pre_IP post_input post_IP  \
-e -p 16 \
--minMappingQuality 20 --skipZeros \
--numberOfSamples 100000 \
-T "Fingerprints of Cut&Run"  \
--plotFile $OUT/hg19.Cut.and.Run.fingerprints.png \
--outRawCounts $OUT/hg19.Cut.and.Run.fingerprints.tab \
--outQualityMetrics $OUT/hg19.Cut.and.Run.QualityMetrics.txt







