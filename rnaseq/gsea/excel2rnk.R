# convert excel to rnk for GSEA
# suitable for DESeq2 output of Rui
# Sort by log2FoldChange_shrinked

excel <- "../DESeq2/K2N_vs_AT.deseq2.xlsx"
outname <- "K2N_vs_AT.rnk"

df <- readxl::read_xlsx( excel, na='-')  # na term important
df <- data.frame(df)  #important
df <- subset(df, select = c("Name","log2FoldChange_shrinked"))
colnames(df) <- c("# Name","log2FoldChange_shrinked")
df <- df[order(df$log2FoldChange_shrinked), ]
head(df)

write.table(df, outname, row.names = F, quote = F, sep='\t')
