SAMPLES=config['SAMPLES']
MODE=config['MODE']


rule fastqc:
    input:
        (expand("fastq/{sample}.{r}.fastq.gz", sample=SAMPLES, r=["R1", "R2"])
            if MODE == 'PE'
            else expand("fastq/{sample}.fastq.gz", sample=SAMPLES))
    output:
        "results/fastqc/multiqc_report.html"
    log:
        "log/fastqc/fastqc.log"
    params:
        odir="results/fastqc"
    threads:
        8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        mkdir -p {params.odir}
        mkdir -p {params.odir}/details
        fastqc -t {threads} {input} -o {params.odir}/details &> {log}
        multiqc {params.odir}/details -o {params.odir} &>> {log}
        """


rule trimmed_reads_fastqc:
    input:
        (expand("results/trimmed_reads/{sample}.{r}.fastq.gz", sample=SAMPLES, r=["R1", "R2"])
            if MODE == 'PE'
            else expand("results/trimmed_reads/{sample}.fastq.gz", sample=SAMPLES))
    output:
        "results/trimmed_reads_fastqc/multiqc_report.html"
    log:
        "log/trimmed_reads_fastqc/fastqc.log"
    params:
        odir="results/trimmed_reads_fastqc"
    threads:
        8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        mkdir -p {params.odir}
        mkdir -p {params.odir}/details
        fastqc -t {threads} {input} -o {params.odir}/details &> {log}
        multiqc {params.odir}/details -o {params.odir} &>> {log}
        """
