#!/bin/bash
#BSUB -J "fastqc.PE[1-4]"
#BSUB -P "Jun"
#BSUB -R rusage[mem=20000]
#BSUB -n 8
#BSUB -q short
#BSUB -W 4:00
##BSUB -N

source config_2018_02_16.sh

cd $workingDir

mkdir -p scripts
mkdir -p scripts/log
mkdir -p scripts/log/fastqc


#BSUB -o /project/umw_nathan_lawson/Jun/2018_02_16_NR2F2_CHIPSEQ/scripts/log/fastqc/PE.out.%J.%I.txt
#BSUB -e /project/umw_nathan_lawson/Jun/2018_02_16_NR2F2_CHIPSEQ/scripts/log/fastqc/PE.err.%J.%I.txt




module load fastqc/0.11.5

# Now let's keep track of some information just in case anything goes wrong

echo "=========================================================="
echo "Starting on       : $(date)"
echo "Running on node   : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID    : $J"
echo "Current job index  : $I"
echo "=========================================================="

mkdir -p $webDir/fastqc_PE
mkdir -p $webDir/fastqc_PE/${sample[$i]}_R1

fastqc $source2/${fq1[$i]} --outdir $webDir/fastqc_PE/${sample[$i]}_R1

mkdir -p $webDir/fastqc_PE/${sample[$i]}_R2

fastqc $source2/${fq2[$i]} --outdir $webDir/fastqc_PE/${sample[$i]}_R2



#mkdir -p fastqc/${sample[$i]}_R2

#fastqc $source/${sample[$i]}_*_R2_*.fastq.gz	--outdir fastqc/${sample[$i]}_R2


echo "=========================================================="
echo "Finished on : $(date)"
echo "=========================================================="



