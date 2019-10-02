# convert excel to rnk for GSEA
# suitable for DESeq2 output of Rui
# Sort by log2FoldChange_shrinked

excel <- "../../../DESeq2_FC1.2/O4_vs_DM4.deseq2.xlsx"
outname <- "O4_vs_DM4.rnk"

df <- readxl::read_xlsx( excel, na='-')  # na term important
df <- data.frame(df)  #important
df <- subset(df, select = c("Name","log2FoldChange_shrinked"))
colnames(df) <- c("# Name","log2FoldChange_shrinked")
df <- df[order(df$log2FoldChange_shrinked), ]
head(df)

write.table(df, outname, row.names = F, quote = F)
