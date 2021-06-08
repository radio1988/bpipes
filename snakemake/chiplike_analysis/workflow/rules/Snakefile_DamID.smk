# June2, 2021
# ChIPSeq like data analysis pipeline
# Functions: 
    # fastqc
    # map, filter, bam_qc
    # call peak: sample level and contrast based (read contras.csv and meta.csv)
    # bigwig (macs2_DamID signal based)
# Env: 
    # source activate chiplike
# Requirements
    # inputs in ./fastq/
    # named as {sample}.{R1,R2}.{fastq,fq}.gz
    # e.g. A.R1.fastq.gz A.R2.fastq.gz B...
    # good format of meta.csv and contrast.csv, matching SAMPLES in config.yaml


from snakemake.utils import min_version
from modules import parse_meta_contrast, get_treatment_bams_from_contrast, \
    get_control_bams_from_contrast, get_contrast_name_from_contrast, \
    get_treat_pileup_bdg_names_from_contrasts, get_treat_pileup_bw_names_from_contrasts, \
    get_control_lambda_bw_names_from_contrasts, get_meme_split_outname_from_contrasts


### parse and prepare
min_version("5.17.0")

configfile: "config/config.yaml"
DATA_TYPE=config['DATA_TYPE']
MODE=config['MODE']
SAMPLES=config["SAMPLES"]
GENOME=config["GENOME"]
INDEX=GENOME+".sa"
MQ_MIN=config["MQ_MIN"]
BIN_SIZE=config["BIN_SIZE"]
GSIZE=config["GSIZE"]
SizeFile=config["SizeFile"]
MEME_DB=config["MEME_DB"]
PEAK_WIDTH=config['PEAK_WIDTH']
CHRS=config['CHRS']
minFragmentLength=config['minFragmentLength']
maxFragmentLength=config['maxFragmentLength']

#o=parse_meta_contrast(fmeta=pwd+"/config/meta.csv", fcontrast=pwd + "/config/contrast.csv") 
# print("parse_meta_contrast_obj:", vars(o))
# {'contrast2contrast_name': {'contrast1': 'G1_vs_ctrl', 'contrast2': 'G2_vs_ctrl', 'contrast3': 'G1_G2_vs_ctrl'}, 
# 'contrast2treatmentSamples': {'contrast1': ['2-1_S2', '2-2_S3', '2-3_S4'], 'contrast2': ['3-1_S5', '3-2_S6', '3-3_S7']}, 
# 'contrast2controlSamples': {'contrast1': ['1-2_S1'], 'contrast2': ['1-2_S1'], 'contrast3': ['1-2_S1']}}

ruleorder: meme_neibour_chr_split > meme_neibour > contrast_treat_pileup_bw > contrast_control_lambda_bw > sample_bdg2bw > split_fa_by_chr > get_summit_neighbour
# Learn: To Avoid AmbiguousRuleException:
# Rules macs2_DamID_contrast_treat_pileup_bw and macs2_DamID_sample_treat_pileup_bw are ambiguous for the file macs2_sample_level_peaks/1-2_S1_treat_pileup.bw.
# Consider starting rule output with a unique prefix, constrain your wildcards, or use the ruleorder directive.
# Wildcards:
#     macs2_DamID_contrast_treat_pileup_bw: sample=1-2_S1
#     macs2_DamID_sample_treat_pileup_bw: sample=1-2_S1
# Expected input files:
#     macs2_DamID_contrast_treat_pileup_bw: macs2_sample_level_peaks/1-2_S1_treat_pileup.bdg
#     macs2_DamID_sample_treat_pileup_bw: macs2_sample_level_peaks/1-2_S1_treat_pileup.bdgExpected output files:
#     macs2_DamID_contrast_treat_pileup_bw: macs2_sample_level_peaks/1-2_S1_treat_pileup.bw
#     macs2_DamID_sample_treat_pileup_bw: macs2_sample_level_peaks/1-2_S1_treat_pileup.bw

shell.prefix("""
            HOME=/home/rl44w/
            bedGraphToBigWig="singularity exec $HOME/singularity/hand_sandbox.simg bedGraphToBigWig"
            """)


### Workflow
# rule targets:
#     input:
#         # 1. everything listed here will be produced by the pipeline
#         # 2. feed {sample}
#         macs2_sample_level_peaks=expand("macs2_sample_level_peaks/{sample}_peaks.narrowPeak", sample=SAMPLES),
#         macs2_DamID_sample_treat_pileup_bw=expand("macs2_sample_level_peaks/{sample}_treat_pileup.bw", sample=SAMPLES),
#         macs2_contrast_level_peaks=get_treat_pileup_bdg_names_from_contrasts(contrasts=o.contrasts, o=o), 
#         macs2_DamID_contrast_treat_pileup_bw=get_treat_pileup_bw_names_from_contrasts(contrasts=o.contrasts, o=o),
#         macs2_DamID_contrast_control_lambda_bw=get_control_lambda_bw_names_from_contrasts(contrasts=o.contrasts, o=o),
# #        meme=get_meme_outname_from_contrasts(contrasts=o.contrasts, PEAK_WIDTH=PEAK_WIDTH, o=o),
#         meme_split=get_meme_split_outname_from_contrasts(contrasts=o.contrasts, PEAK_WIDTH=PEAK_WIDTH, CHRS=CHRS, o=o),
#         # Learn: Good trick to use tagets input to do contrast2contrast_name and more
#         fastqc="fastqc/multiqc_report.html", # not in main workflow, so list here
#         sorted_reads_bam_qc=expand("sorted_reads_bam_qc/stats/{sample}.stats.txt", sample=SAMPLES),
#         multiqc_sorted_reads="sorted_reads_bam_qc/stats/multiqc_report.html",
#         multiqc_DamID_reads="DamID_reads_bam_qc/stats/multiqc_report.html",
#         DamID_reads_bam_qc=expand("DamID_reads_bam_qc/stats/{sample}.stats.txt", sample=SAMPLES),

#         qc1="DamID_reads_bam_qc/fingerprint.pdf",
#         #qc2="DamID_reads_bam_qc/fragment_size.pdf",
#         qc3="DamID_reads_bam_qc/multiBamSummary.heatmap.pdf",
#         qc4="DamID_reads_bam_qc/multiBamSummary.pca.pdf",
#         qc5=expand("DamID_reads_bam_qc/{sample}.insert_size.pdf", sample=SAMPLES),
#         dag="Workflow_DAG.all.svg"


rule fastqc:
    # don't need input, if you agree on not checking them
    # without output, output will not be created
    output:
        "results/fastqc/multiqc_report.html"
    log:
        "log/fastqc/fastqc.log"
    params:
        mem="1000"
    threads:
        8
    conda:
        "../envs/chiplike.yaml"
    shell:
        # {input/output} don't have to be in command
        # have to load module in one block
        """
        module load fastqc/0.11.5
        mkdir -p results/fastqc
        mkdir -p results/fastqc/details
        fastqc -t {threads} fastq/*q.gz -o results/fastqc/details &> {log}
        multiqc results/fastqc/details -o results/fastqc &>> {log}
        """


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


rule bwa_map:
    # 1min/1M reads with 16 cores
    input:
        index=INDEX,
        r1="fastq/{sample}.R1.fastq.gz",
        r2="fastq/{sample}.R2.fastq.gz",
    output:
        temp("results/mapped_reads/{sample}.bam")
    params:
        mem="1500"  # todo auto adjust based on {threads}, for human need 18G+ 
    threads:
        16
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
    # 2M/min
    input:
        "results/mapped_reads/{sample}.bam"
    output:
        "results/sorted_reads/{sample}.bam"
    params:
        mem="1200"
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
        samtools sort -@ {threads} -m 1G {input} -o {output} &>> {log}
        samtools index {output} &>> {log}
        """


rule markDup:
    # same speed as bwa_map, slow
    input:
        "results/sorted_reads/{sample}.bam"
    output:
        #temp(
        bam=temp("results/markDup/{sample}.bam"),
        metrics="results/markDup/{sample}.markDup_metrics.txt",
        bai=temp("results/markDup/{sample}.bam.bai"),
        #)  # make temp and save storage
    log:
        "log/markDup/{sample}.markDup.log"
    benchmark:
        "log/markDup/{sample}.benchmark"
    conda:
        "../envs/chiplike.yaml"
    threads:
        2
    params:
        mem="16000",  # Used Max of 24G before
    shell:
        """
        module load picard/2.17.8
        PICARD=/share/pkg/picard/2.17.8/picard.jar
        
        java -Xmx30g -XX:ParallelGCThreads=2 -jar $PICARD MarkDuplicates \
        I={input} \
        O={output.bam} \
        M={output.metrics} \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true \
        &> {log}

        samtools index {output.bam} &>> {log}
        """


if DATA_TYPE == 'DamID':
    rule DamID_filter:
        input:
            "results/sorted_reads/{sample}.bam"
        output:
            "results/clean_reads/{sample}.bam"
        log:
            "log/clean_reads/{sample}.log"
        benchmark:
            "log/clean_reads/{sample}.benchmark"
        params:
                mem="16000"
        threads:
            1
        conda:
            "../envs/chiplike.yaml"
        shell:
            """
            # need samtools/1.9
            python workflow/scripts/filter_bam.py {input} {GENOME} GATC {output} &> {log}
            """
elif DATA_TYPE == 'ATAC':
    rule ATAC_filter:
        input:
            "results/markDup/{sample}.bam"
        output:
            "results/clean_reads/{sample}.bam"
        log:
            "log/clean_reads/{sample}.log"
        benchmark:
            "log/clean_reads/{sample}.benchmark"
        threads:
            2
        params:
            mem="8000"
        conda:
            "../envs/chiplike.yaml"
        shell:
            """
                echo 'ATACseq mode'
                echo {input}
                echo 'remove mito reads; keep paired reads with MQ>20 and 38-2000nt fragment size only'
                samtools view -h {input} 2>{log}| perl -lane 'print unless ($F[2] eq {chrM} and $_ != /\@/)' 2>>{log}| awk \'{config[filter]}\' 2>>{log}| $samtools sort -m 8G -o {output}  2>> {log}
                cp {input} {output}

            samtools index {output} &> {log}
            """
elif DATA_TYPE == 'ChIP':
    rule ChIP_filter:
        input:
            bam="results/markDup/{sample}.bam",
            bai="results/markDup/{sample}.bam.bai"
        output:
            bam="results/clean_reads/{sample}.bam",
            bai="results/clean_reads/{sample}.bam.bai"
        log:
            "log/clean_reads/{sample}.log"
        benchmark:
            "log/clean_reads/{sample}.benchmark"
        threads:
            2
        params:
            mem="8000"
        conda:
            "../envs/chiplike.yaml"
        shell:
            """
            ln {input.bam} {output.bam} &> {log}
            ln {input.bai} {output.bai} &> {log}
            """
else: 
    sys.exit("DATA_TYPE error, see config.yaml for details")




rule sorted_reads_bam_qc:
    input:
        bam="results/sorted_reads/{sample}.bam"
    output:
        idxstats="results/sorted_reads_bam_qc/idxstats/{sample}.idxstats.txt",
        flagstat="results/sorted_reads_bam_qc/flagstat/{sample}.flagstat.txt",
        stats="results/sorted_reads_bam_qc/stats/{sample}.stats.txt"
    params:
        mem="3000"
    threads:
        1
    log:
        idxstats="log/sorted_reads_bam_qc/idxstats/{sample}.idxstats.log",
        flagstat="log/sorted_reads_bam_qc/flagstat/{sample}.flagstat.log",
        stats="log/sorted_reads_bam_qc/stats/{sample}.stats.log"
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
        stats=expand("results/sorted_reads_bam_qc/stats/{sample}.stats.txt", sample=SAMPLES),
        idxstats=expand("results/sorted_reads_bam_qc/idxstats/{sample}.idxstats.txt", sample=SAMPLES),
        flagstat=expand("results/sorted_reads_bam_qc/flagstat/{sample}.flagstat.txt", sample=SAMPLES),
    output:
        "results/sorted_reads_bam_qc/stats/multiqc_report.html",
        "results/sorted_reads_bam_qc/idxstats/multiqc_report.html",
        "results/sorted_reads_bam_qc/flagstat/multiqc_report.html",
    params:
        mem="3000"
    threads:
        1
    log:
        "log/sorted_reads_bam_qc/multiqc.log"
    conda:
        "../envs/chiplike.yaml"
    shell:
        """
        multiqc -f {input.stats} -o sorted_reads_bam_qc/stats/ &> {log}
        multiqc -f {input.idxstats} -o sorted_reads_bam_qc/idxstats   &>> {log}
        multiqc -f {input.flagstat} -o sorted_reads_bam_qc/flagstat/ &>> {log}
        """

rule clean_reads_bam_qc:
    input:
        bam="results/clean_reads/{sample}.bam"
    output:
        idxstats="results/clean_reads_bam_qc/idxstats/{sample}.idxstats.txt",
        flagstat="results/clean_reads_bam_qc/flagstat/{sample}.flagstat.txt",
        stats="results/clean_reads_bam_qc/stats/{sample}.stats.txt"
    params:
        mem="2000"
    threads:
        1
    log:
        idxstats="log/clean_reads_bam_qc/idxstats/{sample}.idxstats.log",
        flagstat="log/clean_reads_bam_qc/flagstat/{sample}.flagstat.log",
        stats="log/clean_reads_bam_qc/stats/{sample}.stats.log"
    conda:
        "../envs/chiplike.yaml"
    shell:
        """
        samtools idxstats {input.bam} > {output.idxstats} 2> {log.idxstats} 
        samtools flagstat {input.bam} > {output.flagstat} 2> {log.flagstat} 
        samtools stats {input.bam} > {output.stats} 2> {log.stats} 
        """

rule multiQC_clean_reads:
    input:
        stats=expand("results/clean_reads_bam_qc/stats/{sample}.stats.txt", sample=SAMPLES),
        idxstats=expand("results/clean_reads_bam_qc/idxstats/{sample}.idxstats.txt", sample=SAMPLES),
        flagstat=expand("results/clean_reads_bam_qc/flagstat/{sample}.flagstat.txt", sample=SAMPLES),
    output:
        "results/clean_reads_bam_qc/stats/multiqc_report.html",
        "results/clean_reads_bam_qc/idxstats/multiqc_report.html",
        "results/clean_reads_bam_qc/flagstat/multiqc_report.html",
    params:
        mem="3000"
    threads:
        1
    log:
        "log/clean_reads_bam_qc/multiqc.log"
    conda:
        "../envs/chiplike.yaml"
    shell:
        """
        multiqc -f {input.stats} -o DamID_reads_bam_qc/stats/ &> {log}
        multiqc -f {input.idxstats} -o DamID_reads_bam_qc/idxstats   &>> {log}
        multiqc -f {input.flagstat} -o DamID_reads_bam_qc/flagstat/ &>> {log}
        """

rule plotFingerprint:
    input:
        expand("results/clean_reads/{sample}.bam", sample=SAMPLES)
    output:
        plot="results/clean_reads_bam_qc/fingerprint.pdf",
        txt="results/clean_reads_bam_qc/fingerprint.txt",
    params:
        mem="2000"
    threads:
        6
    log:
        "log/clean_reads_bam_qc/fingerprint.log"
    conda:
        "../envs/chiplike.yaml"
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
        plot="results/clean_reads_bam_qc/fragment_size.pdf",
        txt="results/clean_reads_bam_qc/fragment_size.txt"
    params:
        mem="4000"
    threads:
        4
    log:
        "log/clean_reads_bam_qc/fragment_size.log"
    conda:
        "../envs/chiplike.yaml"
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
        "results/clean_reads_bam_qc/multiBamSummary.npz",
    params:
        mem="3500"
    threads:
        8
    log:
        "log/clean_reads_bam_qc/multiBamSummary.log"
    conda:
        "../envs/chiplike.yaml"
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
        "results/clean_reads_bam_qc/multiBamSummary.npz",
    output:
        "results/clean_reads_bam_qc/multiBamSummary.heatmap.pdf"
    params:
        mem="20000"
    threads:
        1
    log:
        "log/clean_reads_bam_qc/plotCorrelation.log"
    conda:
        "../envs/chiplike.yaml"
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
        "results/clean_reads_bam_qc/multiBamSummary.npz",
    output:
        "results/clean_reads_bam_qc/multiBamSummary.pca.pdf"
    params:
        mem="20000"
    threads:
        1
    log:
        "log/clean_reads_bam_qc/plotPCA.log"
    conda:
        "../envs/chiplike.yaml"
    shell:
        """
        plotPCA \
        --corData {input} \
        --plotFile {output} &> {log}
        """

rule CollectInsertSizeMetrics:
    input:
        "results/clean_reads/{sample}.bam"
    output:
        txt="results/clean_reads_bam_qc/{sample}.insert_size.txt",
        pdf="results/clean_reads_bam_qc/{sample}.insert_size.pdf"
    params:
        mem="16000"
    threads:
        1
    # conda: todo
    #     "../envs/chiplike.yaml"
    envmodules:
        "picard/2.17.8"
    shell:
        """
        PICARD=/share/pkg/picard/2.17.8/picard.jar

        java -Xmx15g -jar $PICARD CollectInsertSizeMetrics \
        I={input} \
        O={output.txt} \
        H={output.pdf} &> {log}
        """ 

if DATA_TYPE == 'DamID' and config['MODE'] == 'SITE':
    rule macs2_DamID_sample_SITE:
        input:
            "results/clean_reads/{sample}.bam"
        output:
            "results/macs2_sample_level_peaks/{sample}_peaks.narrowPeak", 
            temp("results/macs2_sample_level_peaks/{sample}_treat_pileup.bdg"),
            temp("results/macs2_sample_level_peaks/{sample}_control_lambda.bdg")
        params:
            mem="8000"    
        threads:
            4
        conda:
            "../envs/macs2.yaml"
        log:
            "log/macs2_sample_level_peaks/{sample}_peaks.narrowPeak.log"
        benchmark:
            "log/macs2_sample_level_peaks/{sample}_peaks.narrowPeak.benchmark"
        shell:
            """
             macs2 callpeak -t {input} \
             -f BAM --nomodel --shift -60 --extsize 100 -g {GSIZE} -q 0.05 --keep-dup all \
             -n {wildcards.sample} --outdir results/macs2_sample_level_peaks -B &> {log}
            """

    rule macs2_DamID_contrast_SITE:
        """
        For each contrast

        MACS2: 
        will concatenate all control bam files and treatment bam files anyway, 
        so no need to collapse tech-reps
        """
        input:
            treatment=lambda wildcards: get_treatment_bams_from_contrast(contrast=wildcards.contrast, o=o),
            control=lambda wildcards: get_control_bams_from_contrast(contrast=wildcards.contrast, o=o),
        output:
            # lambda wildcards: get_contrast_name_from_contrast(contrast=wildcards.contrast)
            "results/macs2_contrast_level_peaks/{contrast}/{contrast_name}_peaks.narrowPeak", # e.g. "macs2_contrast_level_peaks/contrast1/G1_vs_ctrl_peaks.narrowPeak"
            "results/macs2_contrast_level_peaks/{contrast}/{contrast_name}_summits.bed", 
            temp("results/macs2_contrast_level_peaks/{contrast}/{contrast_name}_treat_pileup.bdg"),
            temp("results/macs2_contrast_level_peaks/{contrast}/{contrast_name}_control_lambda.bdg")
        params:
            contrast_name=lambda wildcards: get_contrast_name_from_contrast(contrast=wildcards.contrast, o=o),
            mem="8000"
        threads:
            4        
        conda:
            "../envs/macs2.yaml"
        log:
            "log/macs2_contrast_level_peaks/{contrast}/{contrast_name}.macs2_DamID.log"
        benchmark:
            "log/macs2_contrast_level_peaks/{contrast}/{contrast_name}.macs2_DamID.benchmark"
        shell:
            """
            macs2 callpeak -t {input.treatment} -c {input.control} \
            -f BAM --nomodel --shift -60 --extsize 100 -g {GSIZE} -q 0.05 --keep-dup all \
            -n {params.contrast_name} --outdir results/macs2_contrast_level_peaks/{wildcards.contrast} -B &> {log}
            """
elif (DATA_TYPE == 'DamID' or DATA_TYPE=='ChIP') and config['MODE'] == 'PE':
    rule macs2_DamID_sample_PE:
        input:
            "results/clean_reads/{sample}.bam"
        output:
            "results/macs2_sample_level_peaks/{sample}_peaks.narrowPeak", 
            temp("results/macs2_sample_level_peaks/{sample}_treat_pileup.bdg"),
            temp("results/macs2_sample_level_peaks/{sample}_control_lambda.bdg")
        params:
            mem="8000"    
        threads:
            4
        conda:
            "../envs/macs2.yaml"
        log:
            "log/macs2_sample_level_peaks/{sample}_peaks.narrowPeak.log"
        benchmark:
            "log/macs2_sample_level_peaks/{sample}_peaks.narrowPeak.benchmark"
        shell:
            """
             macs2 callpeak -t {input} \
             -f BAMPE -g {GSIZE} -q 0.05 --keep-dup all \
             -n {wildcards.sample} --outdir results/macs2_sample_level_peaks -B &> {log}
            """

    rule macs2_DamID_contrast_PE:
        """
        For each contrast

        MACS2: 
        will concatenate all control bam files and treatment bam files anyway, 
        so no need to collapse tech-reps
        """
        input:
            treatment=lambda wildcards: get_treatment_bams_from_contrast(contrast=wildcards.contrast, o=o),
            control=lambda wildcards: get_control_bams_from_contrast(contrast=wildcards.contrast, o=o),
        output:
            # lambda wildcards: get_contrast_name_from_contrast(contrast=wildcards.contrast)
            "results/macs2_contrast_level_peaks/{contrast}/{contrast_name}_peaks.narrowPeak", # e.g. "macs2_contrast_level_peaks/contrast1/G1_vs_ctrl_peaks.narrowPeak"
            "results/macs2_contrast_level_peaks/{contrast}/{contrast_name}_summits.bed", 
            temp("results/macs2_contrast_level_peaks/{contrast}/{contrast_name}_treat_pileup.bdg"),
            temp("results/macs2_contrast_level_peaks/{contrast}/{contrast_name}_control_lambda.bdg")
        params:
            contrast_name=lambda wildcards: get_contrast_name_from_contrast(contrast=wildcards.contrast, o=o),
            mem="8000"
        threads:
            4        
        conda:
            "../envs/macs2.yaml"
        log:
            "log/macs2_contrast_level_peaks/{contrast}/{contrast_name}.macs2_DamID.log"
        benchmark:
            "log/macs2_contrast_level_peaks/{contrast}/{contrast_name}.macs2_DamID.benchmark"
        shell:
            """
            macs2 callpeak -t {input.treatment} -c {input.control} \
            -f BAMPE -g {GSIZE} -q 0.05 --keep-dup all \
            -n {params.contrast_name} --outdir results/macs2_contrast_level_peaks/{wildcards.contrast} -B &> {log}
            """
else:
    sys.exit("MODE Error")
    

rule sample_bdg2bw:
    input:
        "results/macs2_sample_level_peaks/{sample}_treat_pileup.bdg"
    output:
        bw="results/macs2_sample_level_peaks/{sample}_treat_pileup.bw",
        sbdg=temp("results/macs2_sample_level_peaks/{sample}_treat_pileup.s.bdg")
    params:
        mem="16000"
    threads:
        1
    log:
        "log/macs2_sample_level_peaks/{sample}_treat_pileup.bw.log"
    benchmark:
        "log/macs2_sample_level_peaks/{sample}_treat_pileup.bw.benchmark"
    priority:
        100
    conda:
        "../envs/chiplike.yaml"
    shell:
        """
        sort -k1,1 -k2,2n {input} > {output.sbdg} 2> {log}
        bedGraphToBigWig="singularity exec $HOME/singularity/hand_sandbox.simg bedGraphToBigWig"
        $bedGraphToBigWig macs2_sample_level_peaks/{wildcards.sample}_treat_pileup.s.bdg {SizeFile} {output.bw} &>> {log}
        """

# rule clean_peaks todo

rule peak2gtf:
    input:
        "results/macs2_sample_level_peaks/{sample}_peaks.narrowPeak"
    output:
        "results/macs2_sample_level_peaks/{sample}_peaks.gtf"
    threads:
        1
    params:
        mem="2000"
    log:
        "log/macs2_sample_level_peaks_peak2gtf/{sample}.log"
    shell:
        """
        perl ..//scripts/peak2gtf.pl {input} > {output} 
        """

rule macs_site_count:
    input:
        bam=expand("results/clean_reads/{sample}.bam", sample=SAMPLES),
        gtf="results/macs2_sample_level_peaks/{sample}_peaks.gtf"
    output:
        "results/macs2_sample_level_peaks/{sample}.count.txt"
    log:
        "log/macs2_sample_level_peaks/{sample}.count.log"
    benchmark:
        "log/macs2_sample_level_peaks/{sample}.count.benchmark"
    threads:
        4
    params:
        mem="8000",
        maxFragmentLength=maxFragmentLength,
        minFragmentLength=minFragmentLength,
        MQ_MIN=MQ_MIN

    conda: 
        '../envs/chiplike.yaml'
    shell:
        """
        featureCounts -a {input.gtf} -o {output} \
        -T {threads} -g gene_id -t peak -s 0 -p -B -C -d {params.minFragmentLength} -D {params.maxFragmentLength} \
        -Q {params.MQ_MIN} --minOverlap 1 --fracOverlap 0 \
        {input.bam} &> {log}
        """

    

rule bamCoverage:
    # for ChIP
    input:
        "results/clean_reads/{sample}.bam"
    output:
        "results/clean_reads_bigWig/{sample}.cpm.bw"
    threads:
        8
    params:
        mem="4000"  # total 6-10G
    log:
        "log/clean_reads_bigWig/{sample}.bamCoverage.log"
    benchmark:
        "log/clean_reads_bigWig/{sample}.bamCoverage.benchmark"
    shell:
        # Aim: same as our downstream filters, extensions
        """
        bamCoverage --bam {input} \
        -o  {output} \
        --numberOfProcessors {threads} \
        --outFileFormat bigwig \
        --normalizeUsing CPM \
        --minFragmentLength {config[minFragmentLength]} \
        --maxFragmentLength {config[maxFragmentLength]} \
        --binSize 10 \
        -e 150 &> {log}
        """

### macs2_contrast_level_peaks ### 

rule get_summit_neighbour:
    input:
        summit="results/macs2_contrast_level_peaks/{contrast}/{contrast_name}_summits.bed",
        genome=GENOME
    output:
        "results/macs2_contrast_level_peaks/{contrast}/{contrast_name}_summits.{width}.fa"
    log:
        "log/macs2_contrast_level_peaks/{contrast}/{contrast_name}_summits.{width}.fa.log"
    threads:
        1
    params:
        mem="6000"
    conda:
        "../envs/chiplike.yaml"
    shell:
        """
        python workflow/scripts/get_summit_neighbour.py {input.genome} {input.summit} {wildcards.width} {output} &> {log}
        """


rule split_fa_by_chr:
    input:
        fasta="results/macs2_contrast_level_peaks/{contrast}/{contrast_name}_summits.{width}.fa"
    output:
        "results/macs2_contrast_level_peaks/{contrast}/by_chr/{contrast_name}_summits.{width}.{chr}.fa"
    log:
        "log/macs2_contrast_level_peaks/{contrast}/by_chr/{contrast_name}_summits.{width}.{chr}.fa.log"
    threads:
        1
    params:
        mem=1000
    conda:
        "../envs/chiplike.yaml"
    shell:
        """
        grep -A 1 '>{wildcards.chr}:' {input} > {output} 2> {log}
        """

rule meme_neibour:
    input: 
        fasta="results/macs2_contrast_level_peaks/{contrast}/{contrast_name}_summits.{width}.fa",
        neg=GENOME,
        db=MEME_DB,
    output: 
        touch("results/macs2_contrast_level_peaks/{contrast}/memechip.{width}/{contrast_name}.finished")
    log:
        "log/macs2_contrast_level_peaks/{contrast}/memechip.{width}/{contrast_name}.log"
    benchmark:
        "log/macs2_contrast_level_peaks/{contrast}/memechip.{width}/{contrast_name}.benchmark"
    params:
        odir="results/macs2_contrast_level_peaks/{contrast}/memechip.{width}/",
        mem="1000",
    threads:
        6
    envmodules:
        "meme/5.0.5"
    shell:
        """
        meme-chip -oc {params.odir} -meme-p {threads} -db {input.db} {input.fasta} &> {log}
        """

rule meme_neibour_chr_split:
    input: 
        fasta=rules.split_fa_by_chr.output,
        neg=GENOME,
        db=MEME_DB,
    output: 
        touch("results/macs2_contrast_level_peaks/{contrast}/memechip_chr.{width}_{chr}/{contrast_name}.finished")
    log:
        "log/macs2_contrast_level_peaks/{contrast}/memechip_chr.{width}_{chr}/{contrast_name}.log"
    benchmark:
        "log/macs2_contrast_level_peaks/{contrast}/memechip_chr.{width}_{chr}/{contrast_name}.benchmark"
    params:
        odir="results/macs2_contrast_level_peaks/{contrast}/memechip_chr.{width}_{chr}/",
        mem="1000",
    threads:
        6
    envmodules:
        "meme/5.0.5"
    shell:
        """
        meme-chip -oc {params.odir} -meme-p {threads} -db {input.db} {input.fasta} &> {log}
        """
        

rule contrast_control_lambda_bw:
    input:
        "results/macs2_contrast_level_peaks/{contrast}/{contrast_name}_control_lambda.bdg"
    output:
        bw="results/macs2_contrast_level_peaks/{contrast}/{contrast_name}_control_lambda.bw",
        sbdg=temp("results/macs2_contrast_level_peaks/{contrast}/{contrast_name}_control_lambda.s.bdg")
    params:
        mem="16000"
    threads:
        1
    priority:
        100
    log:
        "log/macs2_contrast_level_peaks/{contrast}/log/{contrast_name}_control_lambda.bw.log"
    benchmark:
        "log/macs2_contrast_level_peaks/{contrast}/log/{contrast_name}_control_lambda.bw.benchmark"
    conda:
        "../envs/chiplike.yaml"
    shell:
        """
        sort -k1,1 -k2,2n {input} > macs2_contrast_level_peaks/{wildcards.contrast}/{wildcards.contrast_name}_control_lambda.s.bdg 2>> {log}
        bedGraphToBigWig="singularity exec $HOME/singularity/hand_sandbox.simg bedGraphToBigWig"
        $bedGraphToBigWig macs2_contrast_level_peaks/{wildcards.contrast}/{wildcards.contrast_name}_control_lambda.s.bdg \
            {SizeFile} {output.bw} &>> {log}
        """

rule contrast_treat_pileup_bw:
    input:
        "results/macs2_contrast_level_peaks/{contrast}/{contrast_name}_treat_pileup.bdg"
    output:
        bw="results/macs2_contrast_level_peaks/{contrast}/{contrast_name}_treat_pileup.bw",
        sbdg=temp("results/macs2_contrast_level_peaks/{contrast}/{contrast_name}_treat_pileup.s.bdg")
    params:
        mem="16000"
    threads:
        1
    priority:
        100
    log:
        "log/macs2_contrast_level_peaks/{contrast}/{contrast_name}_treat_pileup.bw.log"
    benchmark:
        "log/macs2_contrast_level_peaks/{contrast}/{contrast_name}_treat_pileup.bw.benchmark"
    conda:
        "../envs/chiplike.yaml"
    shell:
        """
        sort -k1,1 -k2,2n {input} > macs2_contrast_level_peaks/{wildcards.contrast}/{wildcards.contrast_name}_treat_pileup.s.bdg 2>> {log}
        bedGraphToBigWig="singularity exec $HOME/singularity/hand_sandbox.simg bedGraphToBigWig"
        $bedGraphToBigWig macs2_contrast_level_peaks/{wildcards.contrast}/{wildcards.contrast_name}_treat_pileup.s.bdg \
            {SizeFile} {output.bw} &>> {log}
        """

rule create_dag:
    params:
        mem="1000"  
        # every job has to have this defined 
        # to use snakemake --cluster 'bsub -q short -R "rusage[mem={params.mem}]" -n {threads}'
    threads:
        1
    output:
        "Workflow_DAG.all.svg"
    log:
        "log/create_dag/Workflow_DAG.all.svg.log"
    shell:
        "snakemake --dag all | dot -Tsvg > {output} 2> {log}"


rule reset:
    shell:
        """
        rm -rf fastqc 
        snakemake --unlock
        """

