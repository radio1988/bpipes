import sys

SAMPLES=config['SAMPLES']
minFragmentLength=config['minFragmentLength']
maxFragmentLength=config['maxFragmentLength']
MQ_MIN=config['MQ_MIN']
BIN_SIZE=config['BIN_SIZE']

def clean_or_sorted_bams_input(wildcards):
   return ["results/"+wildcards.cs_folder+"/"+sample+".bam" for sample in SAMPLES]


def plotFingerprint_params(wildcards):
    if config['MODE'] == "PE":
        return (
            "--minFragmentLength {} --maxFragmentLength {} --extendReads --centerReads".\
            format(config['minFragmentLength'], config['maxFragmentLength'])
        )
    elif config['MODE'] == "SE":
        return " " 
    else:
        sys.exit("MODE Error")

        
rule plotFingerprint_pe:
    input:
        clean_or_sorted_bams_input
    output:
        plot="results/{cs_folder}_qc/fingerprint.pdf",
        txt="results/{cs_folder}_qc/fingerprint.txt",
    params:
        mode=plotFingerprint_params
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2500
    threads:
        6
    log:
        "log/{cs_folder}_qc/fingerprint.log"
    benchmark:
        "log/{cs_folder}_qc/fingerprint.benchmark"
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        plotFingerprint -b {input} \
            --plotFile {output.plot} \
            --outRawCounts {output.txt} \
            --plotTitle "Fingerprint Plot" \
            --smartLabels \
            --binSize {BIN_SIZE} \
            --samFlagInclude 2 \
            -p {threads} {params.mode} &> {log}
        """
        # --samFlagInclude 2: mate properly paired only
        # --extendReads: use mate into

   
rule bamPEFragmentSize:
    input:
        clean_or_sorted_bams_input
    output:
        plot="results/{cs_folder}_qc/fragment_size.pdf",
        txt="results/{cs_folder}_qc/fragment_size.txt"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    threads:
        4
    log:
        "log/{cs_folder}_qc/fragment_size.log"
    benchmark:
        "log/{cs_folder}_qc/fragment_size.benchmark"
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        bamPEFragmentSize \
        -hist {output.plot} \
        --outRawFragmentLengths {output.txt} \
        -T "Fragment Size Distribution" \
        -b {input} \
        --maxFragmentLength 1000 \
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
        --minFragmentLength {minFragmentLength} \
        --maxFragmentLength {maxFragmentLength} \
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


rule insert_size:
    input:
        "results/clean_reads/{sample}.bam"
    output:
        txt="results/clean_reads_qc/insert_size/{sample}.insert_size.txt",
        pdf="results/clean_reads_qc/insert_size/{sample}.insert_size.pdf"
    log:
        'log/clean_reads_qc/insert_size/{sample}.intert_size.log'
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
