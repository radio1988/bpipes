# Notes:
# version 1.0: simple, works
# version 2.0: add args

# read parameters
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop("CMD: Rscript --vanilla rnaseq_normalization_annotation.R ../xx/feature_count.txt ~/dir/gtf.anno_tab.txt out_dir_name\n", call.=FALSE)
} 
count_fname <- args[1]
anno_fname <- args[2]
out_dir <- args[3]
paste("count_table", args[1])
paste("anno_tab", args[2])
dir.create(file.path(out_dir), showWarnings = F)

# load libs
library("scater")
library(WriteXLS)

# read data
df <- read.table(count_fname,
                 sep="\t", header=TRUE, 
                 row.names = 1) # row.name in cts(matrix)

# clean names
colnames(df) <- gsub("\\.bam$", "", colnames(df))
colnames(df) <- gsub("sorted_reads.", "", colnames(df))
row.names(df) <- gsub("\\.*", "", row.names(df))  # for annotation purpose
colnames(df)
dim(df)
cts <- df[, 6:ncol(df)]
cts <- as.matrix(cts)

## Filter Based on Expression
expression_filter <- rowSums(cts) >= 1  # default 10
cts <- cts[expression_filter, ]
df <- df[expression_filter, ]

## read anno
anno <- read.table(
    anno_fname, 
    header=T)
anno$Gene <- gsub("\\.*", "", anno$Gene)  # for annotation purpose
dim(anno)
head(anno)

# TPM calculation
tpm <- calculateTPM(cts, df$Length)
tpm <- data.frame(tpm)
colnames(tpm) <- paste("TPM", colnames(tpm), sep = ":")
head(tpm)
tpm_out <- merge(anno, tpm, by.x=1, by.y=0)
head(tpm_out)
WriteXLS(x = tpm_out, 
         ExcelFileName = paste(out_dir,'TPM.xlsx', sep='/'), row.names = F, SheetNames = 'sheet1', na = '-')

# FPKM calculation
fpkm <- calculateFPKM(cts, df$Length)
fpkm <- data.frame(fpkm)
colnames(fpkm) <- paste("fpkm", colnames(fpkm), sep = ":")
head(fpkm)
fpkm_out <- merge(anno, fpkm, by.x=1, by.y=0)
head(fpkm_out)
WriteXLS(x = fpkm_out, 
         ExcelFileName = paste(out_dir, 'FPKM.xlsx', sep='/'), 
		row.names = F, SheetNames = 'sheet1', na = '-')
