#Generate csv file with readaffy and RMA
#By Lihua Julie Zhu
#on December 7th 2007
#using .CEL file

rm(list=ls())
library(affy)
setwd("~/Documents/ConsultingActivities/MicroarrayExp/MarianWalhout");
require(affy)
Data = ReadAffy(celfile.path="embryo")
#Data<-ReadAffy(filenames=targets$FileName,celfile.path="CEL");
#par(mfrow=c(1,2));

temp = unlist(strsplit(sampleNames(Data), "\\."))

sampleNames(Data) = cbind(temp[1],temp[3], temp[5], temp[7], temp[9], temp[11])

boxplot(Data, col=c(2,3,4,5,6,7))
#boxplot(Data, col=c(1,2,3,4));
eset<-rma(Data);
write.exprs(eset, file="dataRMA.txt");
pcNorm <-read.table("dataRMA.txt", header=TRUE, sep="\t", dec=".");

colnames(pcNorm)[1] = "Probe"

sampleNames(Data)
slotNames(Data)

library("simpleaffy")
Data.qc <- qc(Data)
avbg(Data.qc) #comparable

#Data.qc <- qc.affy(Data,normalised=NULL,tau=0.015,logged=TRUE,
        cdfn=cleancdfname(cdfName(Data)))

#Data.qc <- qc.affy(Data,normalised=NULL,tau=0.015,logged=TRUE, cdfn=cdfName(Data))


#scaling factor
sfs(Data.qc) #comparable

percent.present(Data.qc) # comparable

ratios(Data.qc)[,1:2] # <3 
spikeInProbes(Data.qc)

## Normalization of the data using MAS5.0
eset.mas5 <- mas5(Data)

##################### Get the P/A call info ##############
APInfo <- mas5calls(Data)
#exprs2excel(APInfo, file="Results/dataMas5_PresentCall.csv")

#setwd('./Results')
#exprs2excel(eset.mas5, file="dataMas5.csv")
#write.exprs(eset.mas5, file="Results/dataMas5.csv")

slotNames(APInfo)

present.call <- exprs(APInfo)
colnames(present.call) = paste("PresentCall",colnames(present.call), sep=".");

colnames(present.call)

#boxplot(pcNorm[,2:dim(pcNorm)[2]], col=c(rep(2,4),rep(3,3),rep(4,4),rep(5,4)),range=0);
boxplot(pcNorm[,2:dim(pcNorm)[2]], col=c(2,3,4,5,6,7),range=0);

#par(mfrow=c(2,2));

#image(Data);

#geneIDS <- "need to put a list of ids in"
library(annaffy)

annotation(eset)

#Symbol = aafSymbol(geneIDs, "zebrafish")

##########################################The following is for exploring purpose#############################
deg <- AffyRNAdeg(Data) 
plotAffyRNAdeg(deg,col=c(2,3,4,5,6,7))
#summaryAffyRNAdeg(deg) 
deg$sample.names

legend(5, 10, deg$sample.names, pch = rep(16, 6), col=c(2,3,4,5,6,7))
#legend(0,46,c("LL_Drosophila_2_1","LL_Drosophila_2_2","LL_Drosophila_2_3","LL_Drosophila_2_4"),pch=rep(16,4),col=c(2,3,4,5))

probeNames(Data)[1:10]

dim(pm(Data));
pm(Data)[1:10,]

dim(mm(Data));
dim(intensity(Data)); #intensity for a given probe of the same cdf type across all chips

geneNames(Data)[1:10]

prenorm<- cbind(as.character(probeNames(Data)),pm(Data));
write.table(prenorm, file="preNorm_PM.csv", sep=",");

normPM <- normalize(Data, method="quantiles")
dim(pm(normPM))
qqplot(pm(Data)[,1],pm(Data)[,2]);
qqplot(pm(normPM)[,1], pm(normPM)[,2]);
qqplot(pm(normPM)[,1],pm(normPM)[,3])
qqplot(pm(Data)[,1],pm(Data)[,3])

postnorm <-  cbind(as.character(probeNames(Data)),pm(normPM));
write.table(postPM, file="postNorm_PM.csv", sep=",");

#write.exprs(eset, file="StatRMA.csv");
