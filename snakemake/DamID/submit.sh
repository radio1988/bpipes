#nohup bash submit.sh &
# bsub -W 48:00 -q long 'source activate damid; bash submit.sh &> submit.log '
source activate damid
snakemake -k -p --jobs 999 --latency-wait 20 \
--cluster 'bsub -q long -o lsf.log -R "rusage[mem={params.mem}]" -n {threads} -R span[hosts=1] -W 24:00'

snakemake -j 1 --report report.html
