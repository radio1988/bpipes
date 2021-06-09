
SAMPLES=config['SAMPLES']


rule fastqc_PE:
    input:
        expand("fastq/{sample}.{r}.fastq.gz", sample=SAMPLES, r=["R1", "R2"])
    output:
        "results/fastqc/multiqc_report.html"
    log:
        "log/fastqc/fastqc.log"
    params:
        mem="1000",
        odir="results/fastqc"
    threads:
        8
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