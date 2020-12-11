#nohup bash submit.sh &
# bsub -W 48:00 -q long 'source activate osr; bash submit.sh &> submit.log '

snakemake -k --jobs 999 --latency-wait 600 \
--cluster 'bsub -q long -o lsf.log -R "rusage[mem={params.mem}]" -n {threads} -R span[hosts=1] -W 336:00'

snakemake -j 1 --report report.html

# Notes
# -k: keep going
