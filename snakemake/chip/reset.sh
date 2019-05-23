rm -rf fastqc bam_qc/ mapped_reads/ sorted_reads/ bigWig/ lsf.log *svg *html log nohup.out
rm -rf .snakemake/
snakemake --unlock
