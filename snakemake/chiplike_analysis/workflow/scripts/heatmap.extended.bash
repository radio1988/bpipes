#!/bin/bash

#cd /project/umw_yong-xu_wang/atacseq_peak_caller/liver/ATAC/overlap/heatmap_default
BW=../../genrich_1bp/SRR891269.bw

# Merged Heatmap
computeMatrix scale-regions -a 1000 -b 1000 \
    -S $BW \
    -R ../Genrich.bed ../MACS2.SITE.bed ../MACS2.BAMPE.bed ../HMMRATAC.bed \
    -p 4 \
    -out matrix.full.mat.gz

plotHeatmap -m matrix.full.mat.gz \
    --zMin 0 --zMax 3 \
    --regionsLabel G MS MP H \
    --verbose \
    -out CombinedHeatmap.E.full.pdf &

# Stratified heatmap
computeMatrix scale-regions -a 1000 -b 1000 \
    -S $BW \
    -R ../consensus.bed ../Genrich_uniq.bed ../MACS2.SITE_uniq.bed ../MACS2.BAMPE_uniq.bed ../HMMRATAC_uniq.bed \
    -p 4 \
    -out matrix.unique.mat.gz

plotHeatmap -m matrix.unique.mat.gz \
    --zMin 0 --zMax 3 \
    --regionsLabel C G MS MP H \
    --verbose \
    -out CombinedHeatmap.E.unique5.pdf


# Separate Heatmaps
for BED in ../*bed
do 
name=$(basename $BED)
out=${name/.bed/.tab.gz}
computeMatrix scale-regions -a 1000 -b 1000 -S $BW -R $BED -out $out &
done

wait

for f in *gz
do 
plotHeatmap -m $f -out ${f/.tab.gz/.pdf} --zMin 0 --zMax 3 &
done


# Discarded parameters
## --smartLabels False \
##     --samplesLabel G MS MP H \
## sampleLabel are for bigwig not bed
