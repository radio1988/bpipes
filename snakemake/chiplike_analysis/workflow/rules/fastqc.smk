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
        "../envs/chiplike.yaml"
    shell:
        # {input/output} don't have to be in command
        # have to load module in one block
        """
        mkdir -p results/fastqc
        mkdir -p results/fastqc/details
        fastqc -t {threads} {input} -o {params.odir}/details &> {log}
        multiqc {params.odir}/details -o {params.odir} &>> {log}
        """