## annotate with official GRCz11
perl gtf_cuffcmp_tmap_anno.pl cuffmerge.GRCz11.gtf vs_official.cuffmerge.GRCz11.gtf.tmap Danio_rerio.GRCz11.92.gtf
## annotate leftover with In-house lifted GTF
perl gtf_cuffcmp_tmap_anno.pl cuffmerge.GRCz11.gtf.vs_official.cuffmerge.GRCz11.gtf.tmap.nonanno vs_lifted.cuffmerge.GRCz11.gtf.tmap Danio_rerio.GRCz11.92.gtf
## create merged annotation table for future reference
paste vs_official.cuffmerge.GRCz11.gtf.tmap.txt vs_lifted.cuffmerge.GRCz11.gtf.tmap.txt > both.tmap.txt
## combine 2 steps of annotation
cat cuffmerge.GRCz11.gtf.vs_official.cuffmerge.GRCz11.gtf.tmap.annotated cuffmerge.GRCz11.gtf.vs_official.cuffmerge.GRCz11.gtf.tmap.nonanno.vs_lifted.cuffmerge.GRCz11.gtf.tmap.annotated > cuffmerge_anno.GRCz11.gtf
## sort and index cuffmerge_anno.GRCz11.gtf
bsub < sort.bsub

sleep 360
## remove big temp files
rm -f cuffmerge.GRCz11.gtf.vs_official.cuffmerge.GRCz11.gtf.tmap.annotated

## change to UCSC GTF format(chromosome naming issue only)
perl gtf_format_to_ucsc.pl cuffmerge_anno.GRCz11.gtf chromInfo.txt
mv output.gtf cuffmerge_anno.GRCz11.ucsc.gtf
bsub < sort_ucsc.bsub
rsync -avh cuffmerge_anno.GRCz11.ucsc.gtf* /nl/umw_nathan_lawson/pub/20180628_rnaseq/gtf/
# created hard link in /home/rl44w/umw_nathan_lawson/Annotation/GRCz11/
