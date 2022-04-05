# bsub -W 48:00 -q long -R select[rh=8] 'source activate snakemake; bash workflow/submit.sh &> submit.log '

source activate snakemake

snakemake -k -p --ri \
--notemp \
--use-conda  --conda-prefix "~/anaconda3/envs/mut_analysis/" \
--use-envmodules \
--ri  --restart-times 1 \
--jobs 99  --latency-wait 20 \
--cluster 'bsub -q long -R select[rh=8] -o lsf.log -R "rusage[mem={resources.mem_mb}]" -n {threads} -R span[hosts=1] -W 24:00' \
--cluster-cancel bkill

snakemake -j 1 --report report.html
