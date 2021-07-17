SAMPLES=config['SAMPLES']


rule sorted_reads_qc:
    input:
        bam="results/sorted_reads/{sample}.bam",
        bai="results/sorted_reads/{sample}.bam.bai"
    output:
        idxstats="results/sorted_reads_qc/idxstats/{sample}.idxstats.txt",
        flagstat="results/sorted_reads_qc/flagstat/{sample}.flagstat.txt",
        stats="results/sorted_reads_qc/stats/{sample}.stats.txt"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 3000
    threads:
        1
    log:
        idxstats="log/sorted_reads_qc/idxstats/{sample}.idxstats.log",
        flagstat="log/sorted_reads_qc/flagstat/{sample}.flagstat.log",
        stats="log/sorted_reads_qc/stats/{sample}.stats.log"
    conda:
        "../envs/chiplike.yaml"
    shell:
        """
        samtools idxstats {input.bam} > {output.idxstats} 2> {log.idxstats} 
        samtools flagstat {input.bam} > {output.flagstat} 2> {log.flagstat} 
        samtools stats {input.bam} > {output.stats} 2> {log.stats} 
        """


rule clean_reads_qc:
    input:
        bam="results/clean_reads/{sample}.bam"
    output:
        idxstats="results/clean_reads_qc/idxstats/{sample}.idxstats.txt",
        flagstat="results/clean_reads_qc/flagstat/{sample}.flagstat.txt",
        stats="results/clean_reads_qc/stats/{sample}.stats.txt"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2000
    threads:
        1
    log:
        idxstats="log/clean_reads_qc/idxstats/{sample}.idxstats.log",
        flagstat="log/clean_reads_qc/flagstat/{sample}.flagstat.log",
        stats="log/clean_reads_qc/stats/{sample}.stats.log"
    conda:
        "../envs/chiplike.yaml"
    shell:
        """
        samtools idxstats {input.bam} > {output.idxstats} 2> {log.idxstats} 
        samtools flagstat {input.bam} > {output.flagstat} 2> {log.flagstat} 
        samtools stats {input.bam} > {output.stats} 2> {log.stats} 
        """


rule multiQC_sorted_reads:
    input:
        stats=expand("results/sorted_reads_qc/stats/{sample}.stats.txt", sample=SAMPLES),
        idxstats=expand("results/sorted_reads_qc/idxstats/{sample}.idxstats.txt", sample=SAMPLES),
        flagstat=expand("results/sorted_reads_qc/flagstat/{sample}.flagstat.txt", sample=SAMPLES),
    output:
        "results/sorted_reads_qc/multiqc_report.html",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 3000
    threads:
        1
    log:
        "log/sorted_reads_qc/multiqc.log"
    conda:
        "../envs/chiplike.yaml"
    shell:
        """
        multiqc -f {input.stats} {input.idxstats} {input.flagstat} -o results/sorted_reads_qc/ &> {log}
        """


rule multiQC_clean_reads:
    input:
        stats=expand("results/clean_reads_qc/stats/{sample}.stats.txt", sample=SAMPLES),
        idxstats=expand("results/clean_reads_qc/idxstats/{sample}.idxstats.txt", sample=SAMPLES),
        flagstat=expand("results/clean_reads_qc/flagstat/{sample}.flagstat.txt", sample=SAMPLES),
        markdup=expand("results/markDup/{sample}.markDup_metrics.txt", sample=SAMPLES)
    output:
        "results/clean_reads_qc/multiqc_report.html"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 3000
    threads:
        1
    log:
        "log/clean_reads_qc/multiqc.log"
    conda:
        "../envs/chiplike.yaml"
    shell:
        """
        multiqc -f {input.stats} {input.idxstats} {input.flagstat} {input.markdup} \
        -o results/clean_reads_qc/ &> {log}

        """
