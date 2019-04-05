# Goal: from sample id to run id. 
# Sometimes one sample multiple runs, e.g. SRS454725


# User params
wdir <- '/home/rl44w/mount/dohoon/56samples/download/ena/sample_list'  # sample_list: samples->runs
dbdir <- '/project/umw_dohoon_kim/Rui/56samples/download'
input <- 'GSE48213.56.samples.txt'
output1 <- 'GSE48213.56.samples2runs.txt'
output2 <- 'GSE48213.56.runs.txt'

# process parameters
dbfile <- file.path(dbdir, 'SRAmetadb.sqlite')
input <- file.path(wdir, input)
output1 <- file.path(wdir, output1)
output2 <- file.path(wdir, output2)
#dir.create(odir, showWarnings=F)


# First, set up the library and sqlite connection
library(SRAdb)

if (file.exists(dbfile)){
    sqlfile = dbfile
}else{
    sqlfile = getSRAdbFile(destdir=dbdir)
}

dbcon = dbConnect(dbDriver('SQLite'),sqlfile)


# convert sample ids to run ids
samples <- read.table(input)$V1
runs <- sraConvert(samples, 'run', dbcon)
write.table(runs, output1, quote=F, col.names=F, row.names=F)
write.table(runs$run, output2, quote=F, col.names=F, row.names=F)
