# merge deseq2 and chippeakanno

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

# params
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	chippeakannoFile <- "results/narrow_peaks_contrast_level/c1/HA_vs_IgG_clean.real.narrowPeak.full_anno.xlsx"
  deseq2File <- "results/narrow_peaks_contrast_level/c1/HA_vs_IgG_clean.t_vs_c.xlsx"
  outputFile <- "results/narrow_peaks_contrast_level/c1/HA_vs_IgG_clean.real.narrowPeak.final_anno.xlsx"
  }else{
   chippeakannoFile <- args[1] # ANNO
   deseq2File <- args[2] # TAB
   outputFile <- args[3]
 }




anno <- readExcel(chippeakannoFile)
tab <- readExcel(deseq2File)
merged <- merge(anno, tab, by='peak_id', all.x=T)
merged <- merged[order(merged$pvalue), ]
writeExcel(merged,    outputFile) 
 