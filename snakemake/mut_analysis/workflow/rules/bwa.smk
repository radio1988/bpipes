GENOME=config["GENOME"]
INDEX=GENOME+".sa"
MODE=config['MODE']

def get_bwa_input(wildcards):
    if MODE == 'PE':
        return ["results/trimmed_reads/{sample}.R1.fastq.gz", "results/trimmed_reads/{sample}.R2.fastq.gz"]
    else:
        return "results/trimmed_reads/{sample}.fastq.gz"

rule bwa_index:
    input:
        GENOME
    output:
        INDEX
    log:
        "log/bwa_index.log"
    benchmark:
        'log/bwa_index.benchmark'
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 5000
    threads:
        8
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        bwa index -a bwtsw {input} &> {log}
        """

rule bwa_map:
    '''
    1min/1M reads with 16 cores
    '''
    input:
        index=INDEX,
        reads=get_bwa_input
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
        "../envs/bwa.yaml"
    shell:
        """
        bwa mem -t {threads} {GENOME} {input.reads} 2> {log} | \
        samtools view -Sb -1 -@ 2 - -o {output} &>> {log}
        """



rule bam_sort_index:
# todo: remove temp files, which cause problems when re-run failed submissions
    # 2M/min
    input:
        "results/mapped_reads/{sample}.bam"
    output:
        bam="results/sorted_reads/{sample}.bam",
        bai="results/sorted_reads/{sample}.bam.bai"
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
