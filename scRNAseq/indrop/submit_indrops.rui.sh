logdir="logs/" # path for writing log files
yaml="rui_analysis.yaml"
nlibs=10 # number of libraries to process
nlanes=4


######################################################
# In general, no need to change anything below here. #
######################################################


N1=$nlanes
N3=$nlibs # number of workers for "sort" step
N4=200 # number of workers for "quantify" step


# bug recently introduced that messes up multiple
# workers in the aggregate step
# N5=$nlibs # number of workers for "aggregate" step


module load bcl2fastq/1.8.4 bowtie/1.2.2 jdk/1.8.0_171 RSEM/1.3.0 samtools/1.3
module load gcc/8.1.0


source activate py27
indrops=/project/umw_merav_socolovsky/scRNAseq/tools/indrops/indrops.py
myPython=python


mkdir -p ${logdir}

# filter
bsub -q short -W 4:00 -n 2 -R rusage[mem=3000] -R "span[hosts=1]" -J "filter[1-$N1]" -o ${logdir}filter.%I.out "source ~/scRNAseq.profile; echo \$((\$LSB_JOBINDEX-1)); ${myPython} ${indrops} $yaml filter --total-workers ${N1} --worker-index \$((\$LSB_JOBINDEX-1))"

# abundant_barcode
bsub -w "done(filter)" -q short -W 1:00 -n 1 -R rusage[mem=2000] -R "span[hosts=1]" -J "abundant_barcode" -o ${logdir}abundant_barcode.%I.out "source ~/scRNAseq.profile; ${myPython} $indrops $yaml identify_abundant_barcodes"

# sort
bsub  -w "done(abundant_barcode)" -q long -W 12:00 -n 2 -R rusage[mem=8000] -R "span[hosts=1]" -J "sort[1-$N3]" -o ${logdir}sort.%I.out "source ~/scRNAseq.profile; ${myPython} $indrops $yaml sort --total-workers ${N3} --worker-index \$((\$LSB_JOBINDEX-1))"

# quant
bsub -w "done(sort)" -q short -W 4:00 -n 2 -R rusage[mem=3000] -R "span[hosts=1]" -J "quant[1-$N4]" -o ${logdir}quant.%I.out "source ~/scRNAseq.profile; ${myPython} $indrops $yaml quantify --min-reads 0 --min-counts 0 --total-workers ${N4} --worker-index \$((\$LSB_JOBINDEX-1)) --no-bam"

# aggregate
bsub -w "done(quant)"  -q short -W 2:00  -n 1 -R rusage[mem=8000] -R "span[hosts=1]" -J "aggregate" -o ${logdir}aggregate.%I.out "source ~/scRNAseq.profile; ${myPython} $indrops $yaml aggregate --no-bam"

# trick: \$((\$LSB_JOBINDEX-1)) in LSF
