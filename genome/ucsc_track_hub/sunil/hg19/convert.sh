gffread lincRNA.gff3 -T -o lincRNA.gtf # bioconda cufflinks

rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ ucsc_tools/
 ~/bin/ucsc_tools/gtfToGenePred lincRNA.gtf lincRNA.genePred
 ~/bin/ucsc_tools/genePredToBigGenePred lincRNA.genePred lincRNA.BigGenePred
curl -O https://genome.ucsc.edu/goldenPath/help/examples/bigGenePred.as
~/bin/ucsc_tools/bedToBigBed -type=bed12+8 -tab -as=bigGenePred.as lincRNA.BigGenePred http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes lincRNA.bb

