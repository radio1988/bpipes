GENOME=config["GENOME"]
INDEXFLAG=GENOME+'.bowtie2.done'
MODE=config['MODE']


rule bowtie2_index:
    input:
        GENOME
    output:
        INDEXFLAG
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
        touch {output}
        """

rule bowtie2:
    # 1min/1M reads with 16 cores
    input:
        genome=GENOME,
        index=INDEXFLAG,
        reads=(["fastq/{sample}.R1.fastq.gz", "fastq/{sample}.R2.fastq.gz"]
                  if MODE == 'PE'
                  else "fastq/{sample}.fastq.gz")
    output:
        bam="results/sorted_reads/{sample}.bam",
        bai="results/sorted_reads/{sample}.bam.bai"
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
        bowtie2 -x {input.genome} -p {threads} {params.reads} | \
        samtools sort -@ 2 -m 1G -O BAM -o {output} &>> {log}
        samtools index {output} &>> {log}
        """
