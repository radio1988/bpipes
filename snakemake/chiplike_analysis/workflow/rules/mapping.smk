GENOME=config["GENOME"]
INDEX=GENOME+".sa"


rule bwa_index:
    input:
        GENOME
    output:
        INDEX
    log:
        "log/bwa_index.log"
    params:
        mem="8000"
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
    params:
        mem="3000"  # todo auto adjust based on {threads}, for human need 18G+ 
    threads:
        8
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


rule samtools_sort_index:
# todo: remove temp files, which cause problems when re-run failed submissions
    # 2M/min
    input:
        "results/mapped_reads/{sample}.bam"
    output:
        bam=temp("results/sorted_reads/{sample}.bam"),
        bai=temp("results/sorted_reads/{sample}.bam.bai")
    params:
        mem="2500"
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
        echo sorting finished
        sleep 10
        echo indexing started
        samtools index {output.bam} &>> {log}
        """


# rule bwa_mem_pe:
#     input:
#         reads=["fastq/{sample}.R1.fastq.gz", "fastq/{sample}.R2.fastq.gz"]
#     output:
#         "mapped/{sample}.bam"
#     log:
#         "logs/mapped/bwa.{sample}.log"
#     benchmark:
#         "logs/mapped/bwa.{sample}.benchmark"
#     params:
#         index=config['genome'],
#         extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
#         sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
#         sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
#         sort_extra="-m 2G -@ 2"            # Extra args for samtools/picard.
#     threads: 12
#    # resources:
#    #     mem_mb=lambda wildcards, attempt: attempt * 
#     wrapper:
#         "0.72.0/bio/bwa/mem"

# rule samtools_index:
#     input:
#         "mapped/{sample}.bam"
#     output:
#         "mapped/{sample}.bam.bai"
#     params:
#         "" # optional params string
#     wrapper:
#         "0.72.0/bio/samtools/index"