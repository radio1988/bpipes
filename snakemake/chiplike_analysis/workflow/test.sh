#nohup bash submit.sh &
# bsub -W 48:00 -q long 'source activate damid; bash submit.sh &> submit.log '
source activate snakemake
snakemake -p --notemp --restart-times 0 --ri \
--use-conda  --conda-prefix "~/anaconda3/envs/" \
--use-envmodules -j5
