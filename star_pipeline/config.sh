## genome info
genome_idx='/project/umw_mccb/Rui/genomes/zebra_fish/GRCz11/in_house_merged/star_idx'
genome_fa='/project/umw_mccb/Rui/genomes/zebra_fish/GRCz11/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa'
gtf='/project/umw_mccb/Rui/genomes/zebra_fish/GRCz11/in_house_merged/cuffmerge_anno.GRCz11.gtf'

## aligner info
bam_dir=star

## sample info
fq_dir='../fastq/renamed'

# cd $fq_dir
# ls *R1*
list=('a' 'b' 'c')
R1s=(a_neg.Li1554_R1.fastq.gz  a_pos.Li1555_R1.fastq.gz  neg_neg.Li1558_R1.fastq.gz  pos_pos.Li1560_R1.fastq.gz a_neg.Li1557_R1.fastq.gz  a_pos.Li1556_R1.fastq.gz  neg_neg.Li1559_R1.fastq.gz  pos_pos.Li1561_R1.fastq.gz a_neg.Li1679_R1.fastq.gz  a_pos.Li1678_R1.fastq.gz  neg_neg.Li1687_R1.fastq.gz  pos_pos.Li1845_R1.fastq.gz a_neg.Li1683_R1.fastq.gz  a_pos.Li1844_R1.fastq.gz  neg_neg.Li1848_R1.fastq.gz  pos_pos.Li1847_R1.fastq.gz)
echo 'R1s:' ${#R1s[@]} ${R1s[@]}

# ls *R1*|perl -p -e s/R1/R2/ 
R2s=(a_neg.Li1554_R2.fastq.gz a_neg.Li1557_R2.fastq.gz a_neg.Li1679_R2.fastq.gz a_neg.Li1683_R2.fastq.gz a_pos.Li1555_R2.fastq.gz a_pos.Li1556_R2.fastq.gz a_pos.Li1678_R2.fastq.gz a_pos.Li1844_R2.fastq.gz neg_neg.Li1558_R2.fastq.gz neg_neg.Li1559_R2.fastq.gz neg_neg.Li1687_R2.fastq.gz neg_neg.Li1848_R2.fastq.gz pos_pos.Li1560_R2.fastq.gz pos_pos.Li1561_R2.fastq.gz pos_pos.Li1845_R2.fastq.gz pos_pos.Li1847_R2.fastq.gz)
echo 'R2s:' ${#R2s[@]} ${R2s[@]}

# ls *R1*|perl -p -e s/_R1// | perl -p -e s/.fastq.gz//
names=(a_neg.Li1554 a_neg.Li1557 a_neg.Li1679 a_neg.Li1683 a_pos.Li1555 a_pos.Li1556 a_pos.Li1678 a_pos.Li1844 neg_neg.Li1558 neg_neg.Li1559 neg_neg.Li1687 neg_neg.Li1848 pos_pos.Li1560 pos_pos.Li1561 pos_pos.Li1845 pos_pos.Li1847)
echo 'names:' ${#names[@]} ${names[@]}

## LSF tracker
i=$(($LSB_JOBINDEX- 1))

## tools
featureCounts=/project/umw_mccb/software/subreads/subread-1.5.2-source/bin/featureCounts
bamCoverage=/home/rl44w/.local/bin/bamCoverage
hisat2=/project/umw_mccb/software/hisat2-2.1.0/hisat2
