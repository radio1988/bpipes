# For PullDown vs Input Only
# output Any LFC_RAW > 0
# Use colSum(COUNT) as SizeFactor
# Pvalue not very useful, LFC_RAW good for this step
# todo: Each Treatment must > mean(IgG)?


library(DESeq2)
library(ashr) 
library(WriteXLS)

readExcel <- function(fname){
  df <- readxl::read_xlsx( fname, na='NA', sheet=1)  # na term important
  df <- data.frame(df)  #important
  return (df)
}

writeExcel <- function(df, name){
  WriteXLS(
   df, 
   ExcelFileName=name,
   row.names = FALSE, 
   SheetNames = 'sheet1', 
   na = 'NA')
}

df_filter <- function(df, min_rowsum=10, min_colsum=2e4){
  expression_filter <- rowSums(df[, 6:ncol(df)]) >= min_rowsum
  sample_filter <- colSums(df[, 6:ncol(df)]) >= min_colsum

  if (sum(expression_filter) > 2 & sum(sample_filter) >= 4){
    df <- df[expression_filter, c(rep(TRUE,5), sample_filter )]
    print(paste(
      "Removed genes with less than ", min_rowsum, 
      "reads/fragments across all samples", "\n"))
    print(
      paste("Removed samples with less than ", min_colsum, 
      "reads/fragments across all genes", "\n"))
    print("Data dim after filtering:")
    print(dim(df[, 6:ncol(df)]))
    }else{
      print("Library size too small: \
        too few genes/samples would be left after filtering, so skipped filtering.")
      print("Please interpret results with caution")
    }

    return(df)
}

calculateCPM <- function(counts) {
  # michael's version
  # https://support.bioconductor.org/p/91218/
  x <- counts
  return(t(t(x)*1e6/colSums(x)))
}


clean_name <- function(name){
  # group names in contrast
  name <- gsub(" ", "", name)
  name <- gsub("-", ".", name)
  name <- gsub(";$", "", name)
  return (name)
}

read_csv <- function(fname, sep=','){
  df <- read.csv(fname, comment.char='#', sep=sep)
  df <- df[,colSums(is.na(df))<nrow(df)]  # remove columns with all NA
  df <- data.frame(lapply(df, clean_name))
  return(df)
}

get_contrast <- function(contrast.df, i){
  # get names
  name1 <- clean_name(contrast.df[i,2])
  name2 <- clean_name(contrast.df[i,3])
  name <- paste(name1, name2, sep = "_vs_") # HA.Ab1_vs_IgG
#  name <- paste(contrast.df[i,1], name, sep = '.') # c1.HA.Ab1_vs_IgG

  #if (nchar(name) > 100) {name = contrast.df[i,1]} # c1
  
  resultnames <- gsub("group", "", resultsNames(dds)) #  "group.HA.Ab1" "groupHA.Ab1"  "groupIgG"
  poss <- match(name1, resultnames) # 2
  negs <- match(name2, resultnames) # 3
  contrast <- rep(0, length(resultsNames(dds))) # 0, 0, 0
  contrast[poss] <- 1/length(poss)
  contrast[negs] <- -1/length(negs) # 0, 1, -1; or 0, 1/2, 1/2, -1/2, -1/2 if ; used
  print(data.frame(resNames=resultnames, 
   contrast=contrast))

  return(list(contrast=contrast, name=name))
}

# params
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	countFile <- "results/narrow_peaks_contrast_level/K4me3-KO/K4me3-KO_vs_NoAb-KO_count.txt"
  peakFile <- "results/narrow_peaks_contrast_level/K4me3-KO/K4me3-KO_vs_NoAb-KO_clean.narrowPeak"
  outPeakFile <- "results/narrow_peaks_contrast_level/K4me3-KO/K4me3-KO_vs_NoAb-KO_clean.real.narrowPeak"
	metaFile <- 'config/meta.csv'
	contrastFile <- 'config/contrast.csv'
  contrast_name <- 'K4me3-KO'
  outputExcelFile <- "results/narrow_peaks_contrast_level/K4me3-KO/K4me3-KO_vs_NoAb-KO_clean.t_vs_c.xlsx"
  SizeFactorFile <- "results/clean_reads_qc/stats/reads_mapped.txt"
  }else{
   countFile <- args[1]
   peakFile <- args[2]
   outPeakFile <- args[3]
   metaFile <- args[4]
   contrastFile <- args[5]
   contrast_name <- args[6]
   outputExcelFile <- args[7]
   SizeFactorFile <- args[8]
 }

indfilter <- FALSE
cookscutoff <- TRUE


# Prep
odir <- dirname(countFile)

df 	<- 	read.table(countFile, sep="\t", header=TRUE,  
 comment.char = '#', row.names = 1)
colnames(df) <- gsub("\\.bam$", "", colnames(df))
colnames(df) <- gsub("results.clean_reads.", "", colnames(df))
print(head(df))
print(dim(df))
df[, 6:ncol(df)] <- sapply(df[, 6:ncol(df)], as.integer)

peaks <- read.table(peakFile) # V4 is peak_id
print(head(peaks))

meta <- read_csv(metaFile)
print(meta)

contrast.df <- read_csv(contrastFile)
print(contrast.df)

contrast_name <- clean_name(contrast_name)

# Filter
# df <- df_filter(df, min_rowsum=10, min_colsum=2e4)

pdf(file.path(odir, 'lib_size.pdf'))
boxplot(
  log10(df[, 6:ncol(df)]+1), 
  las=2, 
  main = "library size after filtering"
  )
dev.off()

# COUNT.xlsx, CPM.xlsx
COUNT <- data.frame(df[, 6:ncol(df)])
colnames(COUNT) = paste0(colnames(COUNT), '_COUNT')
print(head(COUNT))
# COUNT.out <- merge(anno, COUNT, by.x=1, by.y=0, all.y=T, sort=F)
# writeExcel(COUNT, file.path(odir, "COUNT.xlsx"))

CPM <- calculateCPM(df[, 6:ncol(df)])
colnames(CPM) = paste0(colnames(CPM), '_CPM')
CPM <- round(CPM,1)
print(head(CPM))
# CPM.out <-  merge(anno, CPM, by.x=1, by.y=0, all.y=T, sort=F)
# writeExcel(CPM,    file.path(odir, "CPM.xlsx"))



# DESeq2
## design
meta <- meta[match(colnames(df[, 6:ncol(df)]), meta$sample), ]
print(meta)

coldata <- data.frame(row.names=colnames(df[, 6:ncol(df)]), 
  sample=factor(meta$sample),
  group=factor(meta$group),
  batch=factor(meta$batch)
  )
print(coldata)

writeLines(capture.output(coldata), 
  file.path(odir,"design.txt"))

if (length(levels(meta$batch)) > 1){
  dds <- DESeqDataSetFromMatrix(
    countData = df[, 6:ncol(df)], 
    colData = coldata, 
    design = ~  0 + group + batch)
  }else{
    dds <- DESeqDataSetFromMatrix(
      countData = df[, 6:ncol(df)], 
      colData = coldata, 
      design = ~  0 + group)  # converted to alph-order
  }

 # !!! unique for ChIPSeq (PullDown vs IgG)
 # If use default, will assume overall binding profile same in PullDown and IgG, 
 # make pos (+LFC) to neg (-LFC)

# Size factor
sf <- read_csv(SizeFactorFile, sep='\t')
sf$count <- as.numeric(as.character(sf$count))
print(sf )
colnames(dds)

sf <- sf[match(colnames(dds),sf$file ), ]
sizeFactors(dds) <- sf$count / min(sf$count)
writeLines(capture.output(sizeFactors(dds)), 
  file.path(odir,"size_factor.txt"))

dds <-DESeq(dds)

NormCount <- counts(dds, normalized=TRUE)
colnames(NormCount) = paste0(colnames(NormCount), '_Norm')



## contrast
i = which(contrast.df$type == contrast_name)
i = i[1]

c = get_contrast(contrast.df, i)
print(c)

if (sum(c$contrast != 0) < 2){
  stop("less than 2 groups in contrast, maybe due to filtering")
}

res.lfc <- lfcShrink(
  dds, 
  contrast = c$contrast, 
  type = 'ashr') 

res.p <- results(
  dds, 
  contrast = c$contrast, 
  independentFilter=indfilter, 
  cooksCutoff=cookscutoff) 

RES <- data.frame(
  peak_id = row.names(res.lfc), 
  log2FC_shrinked = res.lfc$log2FoldChange, 
  log2FC_raw = res.p$log2FoldChange, 
  pvalue = res.p$pvalue,
  FDR = res.p$padj
  )

OUT <- merge(RES, CPM,by.x=1, by.y = 0, all.x=T)
OUT <- merge(OUT, COUNT, by.x=1, by.y = 0, all.x=T)
OUT <- merge(OUT, NormCount,  by.x=1, by.y = 0, all.x=T)

writeExcel(OUT,   outputExcelFile) # "c1/HA_vs_IgG_clean.real.broadPeak"
# "results/narrow_peaks_contrast_level/c1/c1.HA.Ab1_vs_IgG.DESeq2.xlsx"
# subset(OUT, peak_id == 'HA_vs_IgG_peak_2393')
OUT.up <- subset(OUT, log2FC_raw > 0) # for PullDown vs IgG
# OUT.down <- subset(OUT, FDR < MAX_FDR & log2FC_raw < -MIN_LFC)
# OUT.sig <- rbind(OUT.up, OUT.down)
# writeExcel(
#   OUT.up,    
#   file.path(odir, paste0(c$name,".DESeq2.", ".up.xlsx"))) 
peaks.up <- peaks[peaks[,4] %in% OUT.up$peak_id,]
write.table(peaks.up, outPeakFile, sep='\t', row.names=F, col.names=F, quote=F)
peaks.down <- peaks[!peaks[,4] %in% OUT.up$peak_id,]
write.table(peaks.down, paste0(outPeakFile, '.down'), sep='\t', row.names=F, col.names=F, quote=F)


# SessionInfo
writeLines(capture.output(sessionInfo()), 
  file.path(odir,"sessionInfo.txt"))
