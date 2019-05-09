nohup snakemake --jobs 999 --latency-wait 14400 \
--cluster 'bsub -q short -R "rusage[mem={params.mem}]" -n {threads} -R span[hosts=1] -W 4:00' &

snakemake --report report.html

