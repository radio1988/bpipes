#!/bin/bash
# run with:
# nohup bash submit.sh &
# OR: bsub -q long -W 144:00 -R rusage[mem=4000] 'bash submit.sh && echo osc'

module purge
module load singularity/singularity-current > nohup.out   2>&1
source activate osr >> nohup.out  2>&1

snakemake -p -k --jobs 999 \
--use-envmodules \
--use-singularity \
--use-conda  --conda-prefix "/project/umw_mccb/OneStopRNAseq/conda/" \
--latency-wait 300 \
--ri --restart-times 1 \
--cluster 'bsub -q short -o lsf.log -R "rusage[mem={resources.mem_mb}]" -n {threads} -R span[hosts=1] -W 4:00' >> nohup.out  2>&1

snakemake --report report.html > report.log  2>&1
