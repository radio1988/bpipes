SAMPLES=config['SAMPLES']

rule plotFingerprint_PE:
    input:
        expand("results/clean_reads/{sample}.bam", sample=SAMPLES)
    output:
        plot="results/clean_reads_qc/fingerprint.pdf",
        txt="results/clean_reads_qc/fingerprint.txt",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2500
    threads:
        6
    log:
        "log/clean_reads_qc/fingerprint.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        plotFingerprint -b {input} \
            --plotFile {output.plot} \
            --outRawCounts {output.txt} \
            --plotTitle "Fingerprint Plot" \
            --smartLabels \
            --minMappingQuality {MQ_MIN} \
            --binSize {BIN_SIZE} \
            --minFragmentLength {config[minFragmentLength]} \
            --maxFragmentLength {config[maxFragmentLength]} \
            --extendReads \
            --centerReads \
            --samFlagInclude 2 \
            -p {threads} &> {log}
        """
        # --samFlagInclude 2: mate properly paired only
        # --extendReads: use mate into


rule bamPEFragmentSize:
    input:
        expand("results/clean_reads/{sample}.bam", sample=SAMPLES)
    output:
        plot="results/clean_reads_qc/fragment_size.pdf",
        txt="results/clean_reads_qc/fragment_size.txt"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    threads:
        4
    log:
        "log/clean_reads_qc/fragment_size.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        bamPEFragmentSize \
        -hist {output.plot} \
        --outRawFragmentLengths {output.txt} \
        -T "Fragment Size Distribution" \
        --maxFragmentLength 2000 \
        -b {input} \
        -p {threads} &> {log}
        """


rule multiBamSummary:
    input:
        expand("results/clean_reads/{sample}.bam", sample=SAMPLES)
    output:
        "results/clean_reads_qc/multiBamSummary.npz",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 3500
    threads:
        8
    log:
        "log/clean_reads_qc/multiBamSummary.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        multiBamSummary bins \
        -b {input} \
        -o {output} \
        --binSize {BIN_SIZE} \
        --smartLabels \
        -p {threads} \
        --minMappingQuality {MQ_MIN} \
        --minFragmentLength {config[minFragmentLength]} \
        --maxFragmentLength {config[maxFragmentLength]} \
        -e \
        --samFlagInclude 2 &> {log}
        """
        

rule plotCorrelation:
    input:
        "results/clean_reads_qc/multiBamSummary.npz",
    output:
        "results/clean_reads_qc/multiBamSummary.heatmap.pdf"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000
    threads:
        1
    log:
        "log/clean_reads_qc/plotCorrelation.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        plotCorrelation \
        -in {input} \
        --corMethod pearson --skipZeros \
        --whatToPlot heatmap \
        -T 'Pearson Corr Between Bins' \
        --removeOutliers \
        -o {output} &> {log}
        """

rule plotPCA:
    input:
        "results/clean_reads_qc/multiBamSummary.npz",
    output:
        "results/clean_reads_qc/multiBamSummary.pca.pdf"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 20000
    threads:
        1
    log:
        "log/clean_reads_qc/plotPCA.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        plotPCA \
        --corData {input} \
        --plotFile {output} &> {log}
        """

# rule CollectInsertSizeMetrics:
#     input:
#         "results/clean_reads/{sample}.bam"
#     output:
#         txt="results/clean_reads_qc/{sample}.insert_size.txt",
#         pdf="results/clean_reads_qc/{sample}.insert_size.pdf"
#     log:
#         'results/clean_reads_qc/{sample}.intert_size.log'
#     resources:
#         mem_mb=lambda wildcards, attempt: attempt * 1600
#     threads:
#         1
#     # conda: todo
#     #     "../envs/chiplike.yaml"  # also need Rscript
#     envmodules:
#         "picard/2.17.8"
#     shell:
#         """
#         module load picard/2.17.8
#         PICARD=/share/pkg/picard/2.17.8/picard.jar

#         java -Xmx15g -jar $PICARD CollectInsertSizeMetrics \
#         I={input} \
#         O={output.txt} \
#         H={output.pdf} &> {log}
#         """ 


rule insert_size:
    input:
        "results/clean_reads/{sample}.bam"
    output:
        txt="results/clean_reads_qc/insert_size/{sample}.insert_size.txt",
        pdf="results/clean_reads_qc/insert_size/{sample}.insert_size.pdf"
    log:
        'results/clean_reads_qc/{sample}.intert_size.log'
    params:
        # optional parameters (e.g. relax checks as below)
        "VALIDATION_STRINGENCY=LENIENT "
        "METRIC_ACCUMULATION_LEVEL=null "
        "METRIC_ACCUMULATION_LEVEL=SAMPLE"
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=2000
    wrapper:
        "v0.75.0/bio/picard/collectinsertsizemetrics"


# todo: multiqc insert size