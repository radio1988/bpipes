GENOME=config["GENOME"]
INDEX=GENOME+".sa"


rule bwa_index:
    input:
        GENOME
    output:
        INDEX
    log:
        "log/bwa_index.log"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 12000
    threads:
        2
    conda:
        "../envs/chiplike.yaml"
    shell:
        """
        bwa index -a bwtsw {input} &> {log}
        """


rule bwa_map_pe:
    # 1min/1M reads with 16 cores
    input:
        index=INDEX,
        r1="fastq/{sample}.R1.fastq.gz",
        r2="fastq/{sample}.R2.fastq.gz",
    output:
        temp("results/mapped_reads/{sample}.bam")
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2000 # human need 18G
    threads:
        12
    log:
        "log/mapped_reads/{sample}.bam.log"
    benchmark:
        "log/mapped_reads/{sample}.bam.benchmark"
    conda:
        "../envs/chiplike.yaml"
    shell:
        """
        bwa mem -t {threads} {GENOME} \
        {input.r1} {input.r2} \
        2> {log}| samtools view -Sb -1 -@ 2 - -o {output} &>> {log}
        """


rule bam_sort_index:
# todo: remove temp files, which cause problems when re-run failed submissions
    # 2M/min
    input:
        "results/mapped_reads/{sample}.bam"
    output:
        bam=temp("results/sorted_reads/{sample}.bam"),
        bai=temp("results/sorted_reads/{sample}.bam.bai")
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2500
    threads:
        4
    log:
        "log/samtools_sort/{sample}.sort.log"
    benchmark:
        "log/samtools_sort/{sample}.sort.benchmark"
    conda:
        "../envs/chiplike.yaml"
    shell:
        """
        samtools --version &> {log}
        samtools sort -@ {threads} -m 2G {input} -o {output.bam} &>> {log}
        samtools index {output.bam} {output.bai} &>> {log}
        """
