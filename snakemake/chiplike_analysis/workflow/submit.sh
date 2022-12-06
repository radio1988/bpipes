#nohup bash submit.sh &
# bsub -W 48:00 -q long -R select[rh=8] 'source activate snakemake7; bash workflow/submit.sh &> submit.log '

source activate snakemake7


## use conda
snakemake -kp --ri \
--notemp \
--use-conda  --conda-prefix "~/anaconda3/envs/chip" \
--conda-frontend mamba \
--ri  --restart-times 1 \
--jobs 99  --latency-wait 20 \
--cluster 'bsub -q long -R select[rh=8] -o lsf.log -R "rusage[mem={resources.mem_mb}]" -n {threads} -R span[hosts=1] -W 144:00'


## for meme, use module
snakemake -k -p --ri \
--notemp \
--use-envmodules \
--ri  --restart-times 1 \
--jobs 99  --latency-wait 20 \
--cluster 'bsub -q long -R select[rh=6] -o lsf.log -R "rusage[mem={resources.mem_mb}]" -n {threads} -R span[hosts=1] -W 144:00'

snakemake -j 1 --report report.html
