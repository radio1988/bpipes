nohup snakemake --jobs 999 --latency-wait 14400 \
--cluster 'bsub -q long -o lsf.log -R "rusage[mem={params.mem}]" -n {threads} -R span[hosts=1] -W 24:00' &

snakemake --report report.html

