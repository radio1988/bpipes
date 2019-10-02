rm -rf fastqc bam_qc/ mapped_reads/ sorted_reads/ bigWig/ lsf.log *svg *html log nohup.out chip_qc/
rm -rf .snakemake/
snakemake --unlock
