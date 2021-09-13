library("ChIPpeakAnno")
library("GenomicFeatures")
library(rtracklayer)  # for import() function
library(WriteXLS)
require(dplyr)

##todo: add intron/exon in anno
##todo: figure out nearest location vs inside 

# params
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
      peakFile <- "results/narrow_peaks_contrast_level/c1/Mes_vs_Ctrl_clean.narrowPeak"
      gtf_file <- "/project/umw_mccb/genome/Mus_musculus_UCSC_mm10/gencode.vM25.primary_assembly.annotation.gtf" # support GenCode format
      CHIPPEAKANNO_MODE <- 'both'
      BIDING_LEFT <- 3000
      BIDING_RIGHT <- 3000
}else{
      peakFile <- args[1]
      gtf_file <- args[2]
      CHIPPEAKANNO_MODE <- args[3]
      BIDING_LEFT <- strtoi(args[4])
      BIDING_RIGHT <- strtoi(args[5])
}

prefix <- peakFile
odir <- dirname(peakFile)
max_fdr <- 0.05 # for peaks
CHIPPEAKANNO_MODE <- 'both'
bindingRegion <- c(-BIDING_LEFT, BIDING_RIGHT)


# Read gtf
gtf.gr <- import(gtf_file)
gene.gr <- gtf.gr[gtf.gr$type=='gene']
gene.gr.df <- as.data.frame(gene.gr)
transcript.gr <- gtf.gr[gtf.gr$type=='transcript']
TxDb <- makeTxDbFromGRanges(gtf.gr)
TxDb.gr <- toGRanges(TxDb)
TxDb  # todo: save and load


# read and filter peak
peaks.df <- read.table(peakFile, header = F)
dim(peaks.df)
colnames(peaks.df) <- c("seqnames", "start", "end", "peak_id", "score", "strand", "signalValue", "pValue", "qValue" )
peaks.df <- subset(peaks.df, qValue> -log10(max_fdr))
dim(peaks.df)
peaks.gr <- makeGRangesFromDataFrame(peaks.df, keep.extra.columns = T)


# transcript level annotation
gr.anno <- annotatePeakInBatch(myPeakList = peaks.gr, 
                               AnnotationData =  gene.gr, # transcript.gr, gene.gr
                               FeatureLocForDistance="TSS",
                               bindingRegion = bindingRegion,
                               output = CHIPPEAKANNO_MODE, # overlapping, nearestLocation
                               select = 'all',
                               ignore.strand = T)
peak.anno.df <- data.frame(gr.anno)
# select/drop columns (affects function collapse_anno) (to remove too much info in output)
drops <- c('shortestDistance')
peak.anno.df <- peak.anno.df[, !(names(peak.anno.df) %in% drops)]
# merge to anno
peak.anno.df$feature <- as.integer(peak.anno.df$feature)
gene.gr.df$feature <- row.names(gene.gr.df)
meta.df <- gene.gr.df[, c('feature', 'gene_name', 'gene_id', 'gene_type')] # gtf specific
merged <- merge (peak.anno.df, meta.df, by='feature', all.x=T)
peak.anno.df<- merged
# Delete duplicated rows (same peak, same gene) (keep the first)
peak.anno.df$dedup_id <- paste(peak.anno.df$peak_id,peak.anno.df$gene_name)
peak.anno.df<-peak.anno.df %>% distinct(dedup_id, .keep_all = TRUE) # keep the first
peak.anno.df.protein_coding <- subset(peak.anno.df, gene_type=='protein_coding')
WriteXLS(x = peak.anno.df,
         ExcelFileName = paste(prefix, 'anno.WithDup.xlsx', sep = "."),
         row.names = F, SheetNames = 'sheet1', na = '-')  # for user
WriteXLS(x = peak.anno.df.protein_coding,
         ExcelFileName = paste(prefix, 'anno.WithDup.protein_coding.xlsx', sep = "."),
         row.names = F, SheetNames = 'sheet1', na = '-')  # for user
# Collapse by peak_id (double check)
RM <- function(x) gsub(",.*", "", x)
collapse_anno <- function(peak.anno.df, outname){
      # collapse rows by peak_id
      peak.anno.df.collapsed <- peak.anno.df %>%
        group_by(peak_id) %>%
        summarise_each(funs(toString))
      # not collapse seqnames, start, end, width, strand, conc ...
      peak.anno.df.collapsed <- data.frame(apply(peak.anno.df.collapsed[,1:10], 2, RM), 
                                           peak.anno.df.collapsed[,11:14])  # number of columns matters
      print(outname)
      print(dim(peak.anno.df.collapsed))
      WriteXLS(x = peak.anno.df.collapsed,
               ExcelFileName = outname,
               row.names = F, SheetNames = 'sheet1', na = '-')  # for user
      return(peak.anno.df.collapsed)
      }
peak.anno.df.collapsed <- collapse_anno(peak.anno.df, 
      paste(prefix, "anno.collapsed.xlsx", sep = "."))
peak.anno.df.protein_coding.collapsed <- collapse_anno(peak.anno.df.protein_coding, 
      paste(prefix, "anno.protein_coding.collapsed.xlsx", sep = "."))
# Full talbe including unannotated
full.collapsed.df <- merge(peaks.df, peak.anno.df.collapsed[, c(1, 11:ncol(peak.anno.df.collapsed))], 
                           by= "peak_id", all.x=T)
WriteXLS(x = full.collapsed.df,
         ExcelFileName = paste(prefix, "full_anno.xlsx", sep = "."), 
         row.names = F, SheetNames = 'sheet1', na = '-')  # for user



# Anno QC
pdf(file.path(odir, 'anno_stats2.pdf'))
barplot(table(gr.anno$insideFeature))
dev.off()

# TSS QC
pdf(file.path(odir, 'TSS.transcripts.pdf'))
binOverFeature(peaks.gr, annotationData=transcript.gr,
               # select = "nearest",
               # PeakLocForDistance="middle",
               # featureSite="FeatureStart",
               radius=10000, nbins=500, FUN=length, 
               errFun=0,
               ylab="count", 
               main="Distribution of aggregated peak numbers around TSS")
dev.off()

# Gene structure QC
aCR<-assignChromosomeRegion(peaks.gr, 
                            nucleotideLevel=F, # peak centric view
                            precedence=c("Promoters", 
                                         "fiveUTRs", 
                                         "threeUTRs", 
                                         "Exons", 
                                         "immediateDownstream",  
                                         "Introns"), # Rui's precedence
                            proximal.promoter.cutoff	<- 3000,
                            immediate.downstream.cutoff	<- 1000, 
                            TxDb=TxDb)
pdf(file.path(odir, 'anno_stats.pdf'))
op <- par(mar=c(11,4,4,2)) # allows the names.arg below the barplot
barplot(aCR$percentage, las=3)
rm(op)
dev.off()

# SessionInfo
writeLines(capture.output(sessionInfo()), 
      file.path(odir,"sessionInfo.txt"))
