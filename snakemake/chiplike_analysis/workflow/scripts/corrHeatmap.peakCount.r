library(PoiClaClu) # 
library(pheatmap)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  
countFile <- 'HA_vs_IgG_count.txt'
outFile <- 'contrast.samplePoisDist.Heatmap.pdf'
}else{
  countFile <- args[1]
  outFile <- args[2]
}

cat(paste(countFile, outFile))

# read featureCount
df <- read.table(countFile, 
                 sep="\t", header=TRUE, comment.char = '#', 
                 row.names = 1) # row.name in cts(matrix)
colnames(df) <- gsub('\\.bam', '', colnames(df))
colnames(df) <- gsub('results\\.clean_reads\\.', '', colnames(df))
cnt <- df[,c(6:ncol(df))]

# pheatmap
poisd <- PoissonDistance(t(cnt))
samplePoisDistMatrix <- as.matrix( poisd$dd ) 
rownames(samplePoisDistMatrix) <- colnames(cnt)
colnames(samplePoisDistMatrix) <- colnames(cnt)

pheatmap(samplePoisDistMatrix, 
         show_colnames = T, show_rownames = F, 
         angle_col = 315, legend = T,
         filename = outFile)
