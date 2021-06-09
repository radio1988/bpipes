#nohup bash submit.sh &
# bsub -W 48:00 -q long 'source activate snakemake; bash workflow/submit.sh &> submit.log '
source activate snakemake
snakemake -k -p --ri \
--use-conda  --conda-prefix "~/anaconda3/envs/" \
--use-envmodules \
--ri  --restart-times 1 \
--jobs 99  --latency-wait 20 \
--cluster 'bsub -q long -o lsf.log -R "rusage[mem={params.mem}]" -n {threads} -R span[hosts=1] -W 24:00'

snakemake -j 1 --report report.html
