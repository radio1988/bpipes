readExcel <- function(fname){
  df <- readxl::read_xlsx( fname, na='NA', sheet=1)  # na term important
  df <- data.frame(df)  #important
  return (df)
}



# params
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  countFile <- "results/broad_peaks_contrast_level/c1/HA_vs_IgG_count.txt"
  annoFile <- "results/broad_peaks_contrast_level/c1/HA_vs_IgG_clean.broadPeak.full_anno.xlsx"
  metaFile <- 'config/meta.csv'
  contrastFile <- 'config/contrast.csv'
  MAX_FDR <- 0.05
  MIN_LFC <- 0.585
  }else{
   countFile <- args[1]
   annoFile <- args[2]
   metaFile <- args[3]
   contrastFile <- args[4]
   MAX_FDR <- as.numeric(args[5])
   MIN_LFC <- as.numeric(args[6])
 }


 df <- readExcel(annoFile) # peak_id, etc.



  OUT.up <- subset(OUT, FDR < MAX_FDR & log2FC > MIN_LFC)
  OUT.down <- subset(OUT, FDR < MAX_FDR & log2FC < -MIN_LFC)
  OUT.sig <- rbind(OUT.up, OUT.down)
  writeExcel(
    OUT.sig,    
    file.path(odir, paste0(c$name,".DESeq2.", "FDR", MAX_FDR, ".LFC", MIN_LFC, ".sig.xlsx"))) 
