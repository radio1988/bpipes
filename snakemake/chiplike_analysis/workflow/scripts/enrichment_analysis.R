# Read final.xlsx and use gene_name for enrichment analysis, with chippeakanno

# Params
max_FDR = 0.1 # FDR

# params
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  peakAnnoFile <- "results/narrow_peaks_contrast_level/c1/HA_vs_IgG_clean.real.narrowPeak.final_anno.xlsx"
  org_eg_eb <- 'org.Hs.eg.db' # org.Mm.eg.db, org.Hs.eg.db
  } else {
    peakAnnoFile <- args[1]
    org_eg_eb <- args[2]
  }

# packages
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager", repos='http://cran.us.r-project.org')

if (!requireNamespace(org_eg_eb, quietly = TRUE))
BiocManager::install(org_eg_eb, update = FALSE, ask = FALSE)

if (!requireNamespace("KEGG.db", quietly = TRUE))
BiocManager::install("KEGG.db", update = FALSE, ask = FALSE)


library(WriteXLS)
library("ChIPpeakAnno")
library(org_eg_eb, character.only = TRUE)
#library('reactome.db') # fails
library("KEGG.db")

# functions
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

# read
odir <- dirname(peakAnnoFile)
peak.anno.df <- readExcel(peakAnnoFile)
keys <- unique(unlist(strsplit(peak.anno.df$gene_name, ", ")))
write.table(keys, file.path(odir, 'genes.txt'), quote=F, row.names=F, col.names=F)

## GO
print(org_eg_eb)

enriched.GO.both.condense <- getEnrichedGO(
  keys, 
  orgAnn=org_eg_eb, 
  feature_id_type = "gene_symbol", 
  maxP=max_FDR,
  minGOterm=10,
  multiAdjMethod= "BH", 
  condense = T  # all genes in one line
  )

## KEGG
kegg <- getEnrichedPATH(
  keys, 
  orgAnn=org_eg_eb, 
  pathAnn = "KEGG.db",
  feature_id_type = "gene_symbol", 
  maxP=max_FDR,
  minPATHterm=10,
  multiAdjMethod= "BH"
  )

kegg <- condenseMatrixByColnames(as.matrix(kegg), "path.id")

#head(kegg)

## OUTPUT
paste("enriched.GO.both.condense$bp:", dim(enriched.GO.both.condense$bp)[1])
if (dim(enriched.GO.both.condense$bp)[1] > 0) {
  writeExcel(
    enriched.GO.both.condense, 
    file.path(odir, 'enrichment.GO_biological_process.xlsx')
    )
}


paste("enriched.GO.both.condense$cc:", dim(enriched.GO.both.condense$cc)[1])
if (dim(enriched.GO.both.condense$cc)[1] > 0) {
  writeExcel(
    enriched.GO.both.condense$cc,
    file.path(odir, 'enrichment.GO_cellular_component.xlsx')
    )
}

paste("enriched.GO.both.condense$mf:", dim(enriched.GO.both.condense$mf)[1])
if (dim(enriched.GO.both.condense$mf)[1]){
  writeExcel(
    enriched.GO.both.condense$mf,
    file.path(odir, 'enrichment.GO_mulecular_function.xlsx')
    )
}

# WriteXLS(x = enriched.PATH.both.condense,
#         ExcelFileName = 'reactome.db.FDR0.05.FC1.5.xlsx', row.names = F, SheetNames = 'sheet1', na = '-')
paste("kegg:", dim(kegg)[1])
if (dim(kegg)[1] > 0){
  writeExcel(
    kegg,
    file.path(odir, 'enrichment.KEGG.xlsx')
    )
}



