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
    benchmark:
        'log/bowtie2_index.benchmark'
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    threads:
        8
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        bowtie2-build --version &> {log}
        bowtie2-build --threads 8 --verbose  {input} {input} &>> {log}
        """

rule bowtie2:
    # 1min/1M reads with 16 cores
    input:
        genome=GENOME,
        index=INDEX,
        reads=(["fastq/{sample}.R1.fastq.gz", "fastq/{sample}.R2.fastq.gz"]
                  if MODE == 'PE'
                  else "fastq/{sample}.fastq.gz")
    output:
        temp("results/mapped_reads/{sample}.sam")
    params:
        reads="-1 fastq/{sample}.R1.fastq.gz -2 fastq/{sample}.R2.fastq.gz" \
            if MODE == 'PE' else '-U fastq/{sample}.fastq.gz'
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
        bowtie2 -x {input.genome} -p {threads} {params.reads} > {output}
       # | \
       # samtools sort -@ 2 -m 1G -O BAM -o {output} &>> {log}
       # samtools index {output} &>> {log}
        # bowtie2-2.2.5 and samtools-1.7 incompatible in bioconda for libcrypto.so.1.0.0
        """



rule sam_sort_index:
# todo: remove temp files, which cause problems when re-run failed submissions
    # 2M/min
    input:
        "results/mapped_reads/{sample}.sam"
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
        "../envs/samtools.yaml"
    shell:
        """
        samtools --version &> {log}
        samtools sort -@ {threads} -m 2G {input} -o {output.bam} &>> {log}
        samtools index {output.bam} {output.bai} &>> {log}
        """
