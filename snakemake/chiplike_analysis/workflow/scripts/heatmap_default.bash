cd /project/umw_yong-xu_wang/atacseq_peak_caller/GM12878/overlap/heatmap_default

BW=../../genrich_1bp/SRR891269.bw

# Separate Heatmaps
for BED in ../*bed
do 
name=$(basename $BED)
out=${name/.bed/.tab.gz}
computeMatrix reference-point -S $BW -R $BED --referencePoint center -a 1000 -b 1000 -out $out &
done

wait

for f in *gz
do 
plotHeatmap -m $f -out ${f/.tab.gz/.pdf} --zMin 0 --zMax 3 &
done


# Merged Heatmap
computeMatrix reference-point -S $BW \
    -R ../Genrich.bed ../MACS2.SITE.bed ../MACS2.BAMPE.bed ../HMMRATAC.bed \
    --referencePoint center -a 1000 -b 1000 \
    -p 4 \
    -out matrix.full.mat.gz

plotHeatmap -m matrix.full.mat.gz \
    --zMin 0 --zMax 3 \
    --regionsLabel G MS MP H \
    --verbose \
    -out CombinedHeatmap.full.pdf

# Stratified heatmap
computeMatrix reference-point -S $BW \
    -R ../consensus.bed ../Genrich_uniq.bed ../MACS2.SITE_uniq.bed ../MACS2.BAMPE_uniq.bed ../HMMRATAC_uniq.bed \
    --referencePoint center -a 1000 -b 1000 \
    -p 4 \
    -out matrix.unique.mat.gz

plotHeatmap -m matrix.unique.mat.gz \
    --zMin 0 --zMax 3 \
    --regionsLabel C G MS MP H \
    --verbose \
    -out CombinedHeatmap.unique5.pdf

# Discarded parameters
## --smartLabels False \
##     --samplesLabel G MS MP H \
## sampleLabel are for bigwig not bed
