## parameters
genome_idx='/project/umw_mccb/Rui/genomes/zebra_fish/GRCz11/official/star_idx'
genome_fa='/project/umw_mccb/Rui/genomes/zebra_fish/GRCz11/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa'
gtf='/project/umw_mccb/Rui/genomes/zebra_fish/GRCz11/official/Danio_rerio.GRCz11.92.gtf'

bam_dir='star'

## LSF tracker
i=$(($LSB_JOBINDEX- 1))

## tools
featureCounts=/project/umw_mccb/software/subreads/subread-1.5.2-source/bin/featureCounts
bamCoverage=/home/rl44w/.local/bin/bamCoverage
hisat2=/project/umw_mccb/software/hisat2-2.1.0/hisat2
