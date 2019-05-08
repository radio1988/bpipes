nohup snakemake --jobs 999 --latency-wait 14400 \
--cluster 'bsub -q short -o lsf.log -R "rusage[mem={params.mem}]" -n {threads} -R span[hosts=1] -W 4:00' &

