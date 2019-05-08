configfile: "config.yaml"


SAMPLES=config["SAMPLES"]
GENOME=config["GENOME"]
INDEX=config["INDEX"]
GTF=config["GTF"]
STRAND=config["STRAND"]


# load modules (have to use """, to keep in one block)
# - alias does not work, have to use $samstat
shell.prefix("""
            # fastqc
            module load fastqc/0.11.5
            module load python3/3.5.0_packages/multiqc/1.4
            # star
            module load java/1.8.0_77
            module load star/2.5.3a
            #module load samtools/dev-2016_06_06
            samtools="singularity exec $HOME/singularity/hand_sandbox samtools"
            samstat="singularity exec $HOME/singularity/hand_sandbox samstat"
            """)

# Requirements
# inputs in ./fastq/xxx.{fastq,fq}.gz
# named as {sample}.{R1,R2}.{fastq,fq}.gz

# SnakeMake Coding Notes:
# input don't have to be used, just for draw nice DAG

rule all:
    input:
        # 1. everything listed here will be produced by the pipeline
        # 2. feed {sample}
        fastqc="fastqc/multiqc_report.html", # not in main workflow, so list here
        bam_qc=expand("bam_qc/samstat/{sample}.bam.samstat.html", sample=SAMPLES), # feed {samples}
        feature_count=expand("feature_count/counts.gene_id.s{strand}.txt", strand=STRAND), # feed {strand}
        dag="dag.all.svg", # create DAG


rule fastqc:
    # don't need input, if you agree on not checking them
    # without output, output will not be created
    output:
        "fastqc/multiqc_report.html"
    params:
        mem="1000"
    threads:
        6
    log:
        "log/fastqc/fastqc.log"
    shell:
        # {input/output} don't have to be in command
        # have to load module in one block
        """
        mkdir -p fastqc
        mkdir -p fastqc/details
        fastqc -t {threads} fastq/*q.gz -o fastqc/details &> {log}
        multiqc fastqc/details -o fastqc &>> {log}
        """


rule star_map:
    input:
        index=INDEX,
        gtf=GTF,
        r1="fastq/{sample}.R1.fastq.gz",
        r2="fastq/{sample}.R2.fastq.gz",
    output:
        "mapped_reads/{sample}.bam"
    params:
        mem="3000"  # todo auto adjust based on {threads}
    threads:
        12
    log:
        "log/star_map/{sample}.star.log"
    shell:
        # align; rename
        """STAR --runThreadN {threads} \
        --genomeDir {input.index} \
        --sjdbGTFfile {input.gtf} \
        --readFilesCommand zcat \
        --readFilesIn {input.r1} {input.r2} \
        --outFileNamePrefix mapped_reads/{wildcards.sample}. \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 3 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.05 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
        --outSAMstrandField None \
        --outSAMtype BAM Unsorted \
        --quantMode GeneCounts \
        --outReadsUnmapped Fastx \
        &> {log}

        mv mapped_reads/{wildcards.sample}*.out.bam mapped_reads/{wildcards.sample}.bam
        """


rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    params:
        mem="1200"
    threads:
        4
    log:
        "log/samtools_sort/{sample}.sort.log"
    run:
        shell("$samtools sort -@ {threads} -m 1G {input} -o sorted_reads/{wildcards.sample}.bam &> {log}")


rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    params:
        mem="2000"
    threads:
        1
    log:
        "log/samtools_index/{sample}.index.log"
    shell:
        "$samtools index {input} &> {log}"


rule bam_qc:
    input:
        bam="sorted_reads/{sample}.bam",
        bai="sorted_reads/{sample}.bam.bai"
    output:
        "bam_qc/samstat/{sample}.bam.samstat.html"
    params:
        mem="2000"
    threads:
        4
    log:
        idxstats="log/bam_qc/idxstats/{sample}.idxstats.log",
        flagstat="log/bam_qc/flagstat/{sample}.flagstat.log",
        stats="log/bam_qc/stats/{sample}.stats.log",
        samstat="log/bam_qc/samstat/{sample}.samstat.log",

    shell:
        """
        mkdir -p bam_qc 
        mkdir -p bam_qc/idxstats
        mkdir -p bam_qc/flagstat
        mkdir -p bam_qc/stats
        mkdir -p bam_qc/samstat
        $samtools idxstats {input.bam} > bam_qc/idxstats/{wildcards.sample}.idxstats.txt 2> {log.idxstats} &
        $samtools flagstat {input.bam} > bam_qc/flagstat/{wildcards.sample}.flagsat.txt 2> {log.flagstat} &
        $samtools stats {input.bam} > bam_qc/stats/{wildcards.sample}.stats.txt 2> {log.stats} &
        $samstat {input.bam} && mv sorted_reads/{wildcards.sample}*.samstat.html bam_qc/samstat 2> {log.samstat}
        """


rule feature_count:
    input:
        bams=expand("mapped_reads/{sample}.bam", sample=SAMPLES), # for star, faster counting
        gtf=GTF
    output:
        "feature_count/counts.gene_id.s{strand}.txt"
    params:
        mem="4000",
        common="-g gene_id -Q 20 --minOverlap 2 --fracOverlap 0.2",
        pair_end="-p -B -d 50 -D 1000 -C"
    threads:
        4
    log:
        "log/feature_count/counts.gene_id.s{strand}.log"
    shell:
        """
        featureCounts -a {input.gtf} -o {output} \
        -T {threads} \
        {params.common} {params.pair_end} \
        -s {wildcards.strand} \
        {input.bams} &> {log}
        """
        # -p: count Fragments rather than reads for Paired-end reads (remove this for SE data)
        # -C: exclude chimeric (most times, for cancer maybe not include this)
        # -d 50, -D 1000: including PE reads with fragment length in between, which is better than default 50-600
        # -T: num threads
        # -s: strand info, very important; use $i to perform all three possibilities, pick the correct one after counting
        # -Q: min MAPQ, if MAPQ from star, we need to be careful, because star ignores PE information, we might need to add addional step to rescue PE info. (https://github.com/alexdobin/STAR/issues/615)
        # -M: count multiple-mapping reads, based on NH, not useful for RNA-seq, may create misleading summary, by counting multi-mapping reads several times
        # -B: Only count read pairs that have both ends aligned.
        # --fracOverlap 0.2: 20% of read length
        # -–minOverlap 2: 2bp


rule create_dag:
    params:
        mem="1000"  
        # every job has to have this defined 
        # to use snakemake --cluster 'bsub -q short -R "rusage[mem={params.mem}]" -n {threads}'
    threads:
        1
    output:
        "dag.all.svg"
    log:
        "create_dag/dag.all.svg.log"
    shell:
        "snakemake --dag all | dot -Tsvg > {output} 2> {log}"

