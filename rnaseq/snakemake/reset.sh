rm -rf fpkm_tpm/ fastqc/ bigWig/ bam_qc/ mapped_reads/ sorted_reads/ create_dag/ lsf.log *svg *html feature_count log nohup.out
rm -rf .snakemake/
snakemake --unlock
