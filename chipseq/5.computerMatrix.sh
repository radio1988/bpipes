#!/bin/bash
#BSUB -J "5.plot"
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
multiBamSummary=/home/jy91w/.local/bin/multiBamSummary
plotCorrelation=/home/jy91w/.local/bin/plotCorrelation
plotPCA=/home/jy91w/.local/bin/plotPCA
plotCoverage=/home/jy91w/.local/bin/plotCoverage
computeMatrix=/home/jy91w/.local/bin/computeMatrix
plotHeatmap=/home/jy91w/.local/bin/plotHeatmap
plotProfile=/home/jy91w/.local/bin/plotProfile



geneModel=/project/umw_mccb/Jun/genome/genes_bed/ucsc.hg19.knownGenes.bed


ALIGNFOLDER=bam
BWFOLDER=deeptools
OUT=deeptools


workingDir='/project/umw_leslie_shaw/Jun_Allele/Mag/20170720-CUTandRUN'
cd $workingDir

mkdir -p scripts
mkdir -p scripts/log
mkdir -p scripts/log/dpt


#BSUB -o /project/umw_leslie_shaw/Jun_Allele/Mag/20170720-CUTandRUN/scripts/log/pdt/5a.out.%J.%I.txt
#BSUB -e /project/umw_leslie_shaw/Jun_Allele/Mag/20170720-CUTandRUN/scripts/log/dpt/5a.err.%J.%I.txt


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
$computeMatrix scale-regions --numberOfProcessors 16 \
-R  $geneModel  \
-S ${BWFOLDER}/pre_input.bw ${BWFOLDER}/pre_IP.bw ${BWFOLDER}/post_input.bw ${BWFOLDER}/post_IP.bw  \
--samplesLabel pre_input pre_IP post_input post_IP \
-b 5000 -a 5000 --regionBodyLength 5000 \
--skipZeros -o $OUT/hg19_matrix_scaled.gz \
--outFileNameMatrix $OUT/hg19_matrix_scaled.tab \
--outFileSortedRegions $OUT/hg19_matrix_scaled_genes.bed

$plotHeatmap -m $OUT/hg19_matrix_scaled.gz \
--samplesLabel pre_input pre_IP post_input post_IP \
-out $OUT/hg19_heatmap.pdf \
--dpi 100 \
--boxAroundHeatmaps no \
--plotTitle "Cut&Run heatmap"

$plotProfile -m $OUT/hg19_matrix_scaled.gz \
-out $OUT/hg19_profile.pdf \
--samplesLabel pre_input pre_IP post_input post_IP \
--dpi 100 --perGroup --plotTitle "Cut&Run heatmap"

$plotProfile -m $OUT/hg19_matrix_scaled.gz \
-out $OUT/hg19_profile.kmeans.2.pdf \
--samplesLabel pre_input pre_IP post_input post_IP \
--dpi 100 --perGroup --numPlotsPerRow 2 --kmeans 2 \
--plotTitle "Cut&Run heatmap kmeans=2"


#$plotFingerprint \
#-b ${ALIGNFOLDER}/Artery-5.sort.markDup.bam ${ALIGNFOLDER}/Artery-10.sort.markDup.bam ${ALIGNFOLDER}/Artery-15.sort.markDup.bam ${ALIGNFOLDER}/Vein-5.sort.markDup.bam ${ALIGNFOLDER}/Vein-10.sort.markDup.bam ${ALIGNFOLDER}/Vein-15.sort.markDup.bam  \
#--labels Artery-5 Artery-10 Artery-15 Vein-5 Vein-10 Vein-15  \
#-e -p 16 \
#--minMappingQuality 20 --skipZeros \
#--numberOfSamples 100000 \
#-T "Fingerprints of Cut&Run"  \
#--plotFile $OUT/hg19.Cut.and.Run.fingerprints.png \
#--outRawCounts $OUT/hg19.Cut.and.Run.fingerprints.tab \
#--outQualityMetrics $OUT/hg19.Cut.and.Run.QualityMetrics.txt


$multiBamSummary bins --bamfiles ${ALIGNFOLDER}/pre_input.sort.markDup.bam ${ALIGNFOLDER}/pre_IP.sort.markDup.bam ${ALIGNFOLDER}/post_input.sort.markDup.bam ${ALIGNFOLDER}/post_IP.sort.markDup.bam  \
--labels pre_input pre_IP post_input post_IP   \
-bs 10 -p 16 -e \
--minMappingQuality 20  \
-out $OUT/hg19.Cut.and.Run.multiBamSummary.npz \
--outRawCounts $OUT/hg19.Cut.and.Run.multiBamSummary.txt
#
#
$plotCoverage \
--bamfiles ${ALIGNFOLDER}/pre_input.sort.markDup.bam ${ALIGNFOLDER}/pre_IP.sort.markDup.bam ${ALIGNFOLDER}/post_input.sort.markDup.bam ${ALIGNFOLDER}/post_IP.sort.markDup.bam \
--labels pre_input pre_IP post_input post_IP  \
--plotFile $OUT/hg19.Cut.and.Run_coverage.pdf \
-p 16 -e \
-n 1000000 \
--plotTitle "coverage of cut&run" \
--outRawCounts $OUT/hg19.Cut.and.Run_coverage.txt \
--ignoreDuplicates \
--minMappingQuality 20
#
#
#
$plotCorrelation \
-in $OUT/hg19.Cut.and.Run.multiBamSummary.npz \
--labels pre_input pre_IP post_input post_IP  \
--corMethod pearson --skipZeros \
--removeOutliers \
--plotTitle "Pearson Correlation" \
--whatToPlot scatterplot \
-o $OUT/hg19.Cut.and.Run.scatterplot_PearsonCorr.pdf   \
--outFileCorMatrix $OUT/hg19.Cut.and.Run.scatterplot_PearsonCorr.txt


$plotCorrelation -in $OUT/hg19.Cut.and.Run.multiBamSummary.npz \
--labels pre_input pre_IP post_input post_IP  \
--corMethod spearman --skipZeros \
--plotTitle "Spearman Correlation" \
--removeOutliers \
--whatToPlot scatterplot \
-o $OUT/hg19.Cut.and.Run.heatmap_SpearmanCorr.pdf   \
--outFileCorMatrix $OUT/hg19.Cut.and.Run.scatter_SpearmanCorr.txt
#
$plotPCA -in $OUT/hg19.Cut.and.Run.multiBamSummary.npz  \
--labels pre_input pre_IP post_input post_IP  \
--plotTitle "PCA of CUT&RUN"  \
--outFileNameData $OUT/hg19.Cut.and.Run.PCA.data.txt \
-o $OUT/hg19.Cut.and.Run.PCA.pdf





