# Use peak_id in real.narrowPeak to filter summit.bed

# params
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  summitFile <- "results/narrow_peaks_contrast_level/c1/HA_vs_IgG_summits.bed"
  realPeakFile <- "results/narrow_peaks_contrast_level/c1/HA_vs_IgG_clean.real.narrowPeak"
  realSummitFile <-  "results/narrow_peaks_contrast_level/c1/HA_vs_IgG_summits.real.bed"
  }else{
   summitFile <- args[1]
   realPeakFile <- args[2]
   realSummitFile <- args[3]
 }

# Prep
odir <- dirname(summitFile)
realPeaks <- read.table(realPeakFile) # V4 is peak_id
summits <- read.table(summitFile) # V4 peak_id
realSummits <- summits[summits[,4] %in% realPeaks[,4],]
write.table(realSummits, realSummitFile, sep='\t', row.names=F, col.names=F, quote=F)
