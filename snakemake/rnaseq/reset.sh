rm -rf fastqc bam_qc/ mapped_reads/ sorted_reads/ *svg *html feature_count log nohup.out
rm -rf .snakemake/
snakemake --unlock
