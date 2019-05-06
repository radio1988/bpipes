#!/bin/sh
# properties = {"wildcards": ["S2"], "local": false, "rule": "star_map", "threads": 12, "input": ["fastq/S2.fastq.gz", "/project/umw_mccb/genome/Homo_sapiens/ucsc_hg38_primary/gencode.v29.primary_assembly.annotation.fixed.gtf", "/project/umw_mccb/genome/Homo_sapiens/ucsc_hg38_primary/star_idx"], "cluster": {}, "log": [], "output": ["mapped_reads/S2.bam"], "params": {}, "resources": {}, "jobid": 9}
cd /home/rl44w/github/bpipes/snakemake/rnaseq && \
/home/rl44w/miniconda3/envs/py35/bin/python -m snakemake mapped_reads/S2.bam --snakefile /home/rl44w/github/bpipes/snakemake/rnaseq/Snakefile \
--force -j --keep-target-files --keep-shadow --keep-remote \
--wait-for-files /home/rl44w/github/bpipes/snakemake/rnaseq/.snakemake/tmp.ymji6mox fastq/S2.fastq.gz /project/umw_mccb/genome/Homo_sapiens/ucsc_hg38_primary/gencode.v29.primary_assembly.annotation.fixed.gtf /project/umw_mccb/genome/Homo_sapiens/ucsc_hg38_primary/star_idx --latency-wait 5 \
--benchmark-repeats 1 \
\
--force-use-threads --wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
  -p --nocolor \
--notemp --quiet --no-hooks --nolock --printshellcmds  --force-use-threads  --allowed-rules star_map  && touch "/home/rl44w/github/bpipes/snakemake/rnaseq/.snakemake/tmp.ymji6mox/9.jobfinished" || (touch "/home/rl44w/github/bpipes/snakemake/rnaseq/.snakemake/tmp.ymji6mox/9.jobfailed"; exit 1)

