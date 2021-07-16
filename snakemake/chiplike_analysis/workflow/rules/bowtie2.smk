GENOME=config["GENOME"]
INDEX=GENOME+".rev.2.bt2"
MODE=config['MODE']


rule bowtie2_index:
    input:
        GENOME
    output:
        INDEX
    log:
        "log/bowtie2_index.log"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 6000
    threads:
        4
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        bowtie2-build --version > {log}
        bowtie2-build --threads {threads} {input} {input} &>> {log}
        """

rule bowtie2:
    # 1min/1M reads with 16 cores
    input:
        index=INDEX,
        reads=(expand("fastq/{sample}.{r}.fastq.gz", r=["R1", "R2"])
                  if MODE == 'PE'
                  else "fastq/{sample}.fastq.gz")
    output:
        temp("results/mapped_reads/{sample}.bam")
    params:
        reads="-1 fastq/{sample}.R1.fastq.gz -2 fastq/{sample}.R2.fastq.gz" \
            if config['PAIR_END'] else '-U fastq/{sample}.fastq.gz'
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000 # human need 18G
    threads:
        12
    log:
        "log/mapped_reads/{sample}.bowtie2.log"
    benchmark:
        "log/mapped_reads/{sample}.bowtie2.benchmark"
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        bowtie2 --version &> {log}
        bowtie2 -x {input.genome} -p {threads} {params.reads} \
        | samtools sort -@ 2 -m 1G -O BAM -o {output.bam} &>> {log}
        samtools index {output.bam} &>> {log}
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
