args = commandArgs(trailingOnly = T)
 
bam_folder= args[1]
outDir <- args[2]
dir.create(outDir)
sampleName = args[3]
 
treatBAM =  args[4]
inputBAM = args[5]
## window size for tile of genome. Here we set big window size because that
## 1. peaks of NADs are wide;
## 2. data could be smoothed well.
 
ws = as.numeric(args[6]) #50000
step = as.numeric(args[7])
backgroundPercentage =  as.numeric(args[8])
inputBamfile <- file.path(bam_folder, inputBAM)
treatBamfile <- file.path(bam_folder, treatBAM)
 
 
library(GenomicAlignments)
library(GenomeInfoDb)
library(Rsamtools)
#library(rbamtools)
library(csaw)
library(SummarizedExperiment)
#library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Hsapiens.UCSC.hg38)
 
####### use tileCount function
library(NADfinder)
reads = c(inputBamfile, treatBamfile)
all(grepl(".bam$", reads))
 
all(file.exists(paste0(reads, ".bai")))
se <- tileCount(reads=reads,
                genome=Hsapiens,
                windowSize=ws, 
                step=step,
                mode = IntersectionNotStrict,
                dataOverSamples = FALSE)
 
saveRDS(se, file.path(outDir, paste0("1.tileCounts.bsub.", sampleName, ".RDS")))
