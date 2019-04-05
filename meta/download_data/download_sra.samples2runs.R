# Goal: from sample id to run id. 
# Sometimes one sample multiple runs, e.g. SRS454725


# First, set up the library and sqlite connection
library(SRAdb)
sqlfile = getSRAdbFile()
dbcon = dbConnect(dbDriver('SQLite'),sqlfile)

# convert sample ids to run ids
input <- '56samples.sra.txt'
output1 <- '56samples.sra.samples2runs.txt'
output2 <- '56samples.sra.runs.txt'

samples <- read.table(input)$V1
runs <- sraConvert(samples, 'run', dbcon)
write.table(runs, output1, quote=F, col.names=F, row.names=F)
write.table(runs$run, output2, quote=F, col.names=F, row.names=F)
