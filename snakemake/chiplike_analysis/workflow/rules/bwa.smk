rule bwa_mem_pe:
    input:
        reads=["fastq/{sample}.R1.fastq.gz", "fastq/{sample}.R2.fastq.gz"]
    output:
        "mapped/{sample}.bam"
    log:
        "logs/mapped/bwa.{sample}.log"
    benchmark:
        "logs/mapped/bwa.{sample}.benchmark"
    params:
        index=config['genome'],
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra="-m 2G -@ 2"            # Extra args for samtools/picard.
    threads: 12
   # resources:
   #     mem_mb=lambda wildcards, attempt: attempt * 
    wrapper:
        "0.72.0/bio/bwa/mem"

rule samtools_index:
    input:
        "mapped/{sample}.bam"
    output:
        "mapped/{sample}.bam.bai"
    params:
        "" # optional params string
    wrapper:
        "0.72.0/bio/samtools/index"
