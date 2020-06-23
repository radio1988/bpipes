nohup snakemake -k --jobs 999 --latency-wait 1209600 \
--cluster 'bsub -q long -o lsf.log -R "rusage[mem={params.mem}]" -n {threads} -R span[hosts=1] -W 336:00' && snakemake --report report.html &

nohup snakemake --use-conda -k --jobs 999 --latency-wait 14400 \
--cluster 'bsub -q short -o lsf.log -R "rusage[mem={params.mem}]" -n {threads} -R span[hosts=1] -W 4:00' && snakemake --report report.html &

# Notes
# -k: keep going
