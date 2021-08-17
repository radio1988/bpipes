bsub -q long -W 336:00 -n 2 -R span[hosts=1] -R rusage[mem=64000] '
module load blast/2.7.1+; blastn -task blastn -outfmt 7 -word_size 5 -evalue 100000000000000000000 -num_threads 2 -db /home/rl44w/mccb/genome/Homo_sapiens/ucsc_hg38_primary/hg38.primary.fa -query gli.fa > gli.blast7' # max 48 G RAM USED
