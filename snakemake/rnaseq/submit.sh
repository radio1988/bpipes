nohup snakemake --jobs 999 --cluster 'bsub -q short -R "rusage[mem=3000]" -n 12 -R span[hosts=1] -W 4:00' &

