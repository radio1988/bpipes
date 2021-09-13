if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("org.Mm.eg.db", quietly = TRUE))
  BiocManager::install("org.Mm.eg.db")

source("https://bioconductor.org/biocLite.R")
biocLite("KEGG.db")


#library("org.Hs.eg.db")
library("org.Mm.eg.db")
#library('reactome.db')
library("KEGG.db")
library(readxl)
library(WriteXLS)
library("ChIPpeakAnno")


org_eg_eb <- "org.Mm.eg.db"
max_BH = 0.1

peak.anno.df <- read_excel('contrast1.narrow.clean.anno.WithDup.xlsx')
keys <- unique(peak.anno.df$gene_name)

## GO
paste(org_eg_eb)

enriched.GO.both.condense <- getEnrichedGO(
  keys, 
  orgAnn=org_eg_eb, 
  feature_id_type = "gene_symbol", 
  maxP=max_BH,
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
  maxP=max_BH,
  minPATHterm=10,
  multiAdjMethod= "BH"
)

kegg <- condenseMatrixByColnames(as.matrix(kegg), "path.id")

#head(kegg)

## OUTPUT
paste("enriched.GO.both.condense$bp", dim(enriched.GO.both.condense$bp)[1])
WriteXLS(x = enriched.GO.both.condense$bp,
         ExcelFileName = paste('GO_biological_process', max_BH, 'xlsx', sep = '.'), 
         row.names = F, SheetNames = 'sheet1', na = '-')

paste("enriched.GO.both.condense$cc", dim(enriched.GO.both.condense$cc)[1])
WriteXLS(x = enriched.GO.both.condense$cc,
         ExcelFileName = paste('GO_cellular_component', max_BH, 'xlsx', sep = '.'),
         row.names = F, SheetNames = 'sheet1', na = '-')

paste("enriched.GO.both.condense$mf", dim(enriched.GO.both.condense$mf)[1])
WriteXLS(x = enriched.GO.both.condense$mf,
         ExcelFileName = paste('GO_mulecular_function', max_BH, 'xlsx', sep = '.'),
         row.names = F, SheetNames = 'sheet1', na = '-')

# WriteXLS(x = enriched.PATH.both.condense,
#         ExcelFileName = 'reactome.db.FDR0.05.FC1.5.xlsx', row.names = F, SheetNames = 'sheet1', na = '-')

paste("kegg", dim(kegg)[1])
WriteXLS(x = kegg,
         ExcelFileName = paste('KEGG', max_BH, 'xlsx', sep = '.'),
         row.names = F, SheetNames = 'sheet1', na = '-')