#!/bin/bash
#BSUB -J "macs2[1-2]"
#BSUB -R rusage[mem=20000]
#BSUB -n 16
#BSUB -q long
#BSUB -W 24:00
##BSUB -N
#BSUB -P "Nathan"
#BSUB -w 'ended("bowtie")'

#set -u
#set -e
#set -o pipefail
#### -g 1.35e9 for ZV9
###### chromeSize="/project/umw_mccb/genome/danRer7.chrom.sizes"

PICARD=/project/umw_mccb/software/picard2.10.0/picard.jar
macs2=/home/jy91w/.local/bin/macs2
workingDir='/project/umw_leslie_shaw/Jun_Allele/Mag/20170720-CUTandRUN'
cd $workingDir

mkdir -p scripts
mkdir -p scripts/log
mkdir -p scripts/log/macs2


#BSUB -o /project/umw_leslie_shaw/Jun_Allele/Mag/20170720-CUTandRUN/scripts/log/macs2/3A.out.%J.%I.txt
#BSUB -e /project/umw_leslie_shaw/Jun_Allele/Mag/20170720-CUTandRUN/scripts/log/macs2/3A.err.%J.%I.txt


module load python/2.7.9
# module load python/2.7.9_packages/macs2/2.1.0
module load bedtools/2.25.0
module load ucsc_cmdline_util/12_16_2013
chromeSize="/project/umw_mccb/genome/hg19.chrome.size"

i=$(($LSB_JOBINDEX - 1))

sample=(pre_input pre_IP post_input post_IP)
IP=(pre_IP post_IP)
Input=(pre_input post_input)
Signal=(pre_EMT post_EMT)
species=hg19
chromeSize="/project/umw_mccb/genome/hg19.chrome.size"
bdg2bw=/project/umw_mccb/software/bdg2bw
MOUT=macs2
ALIGNFOLDER=bam

module load python/2.7.9

mkdir -p $MOUT
#
$macs2 callpeak -t  ${ALIGNFOLDER}/${IP[$i]}.sort.markDup.bam -c  ${ALIGNFOLDER}/${Input[$i]}.sort.markDup.bam \
-f BAMPE  -g hs -n ${Signal[$i]} --bdg -q 0.05 --call-summits --SPMR --outdir $MOUT
#
$macs2 callpeak --broad -t  ${ALIGNFOLDER}/${IP[$i]}.sort.markDup.bam -c  ${ALIGNFOLDER}/${Input[$i]}.sort.markDup.bam \
-f BAMPE  -g hs -n ${Signal[$i]}.broad -q 0.05 --outdir $MOUT

## Run MACS2 bdgcmp to generate fold-enrichment and logLR track
##  * -m FE means to calculate fold enrichment. Other options can be logLR for log likelihood, subtract for subtracting noise from treatment sample.
## * -p sets pseudocount. This number will be added to 'pileup per million reads' value. You don't need it while generating fold enrichment track because control lambda will always >0. But in order to avoid log(0) while calculating log likelihood, we'd add pseudocount. Because I set precision as 5 decimals, here I use 0.00001.

cd $MOUT
$macs2 bdgcmp -t ${Signal[$i]}_treat_pileup.bdg -c ${Signal[$i]}_control_lambda.bdg -o ${Signal[$i]}_FE.bdg -m FE
$macs2 bdgcmp -t ${Signal[$i]}_treat_pileup.bdg -c ${Signal[$i]}_control_lambda.bdg -o ${Signal[$i]}_logLR.bdg -m logLR -p 0.00001


### Fix the bdg, convert to bigwig

$bdg2bw ${Signal[$i]}_FE.bdg $chromeSize
$bdg2bw ${Signal[$i]}_logLR.bdg $chromeSize








