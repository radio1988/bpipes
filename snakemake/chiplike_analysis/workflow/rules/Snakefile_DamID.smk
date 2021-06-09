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
#         sorted_reads_qc=expand("sorted_reads_qc/stats/{sample}.stats.txt", sample=SAMPLES),
#         multiqc_sorted_reads="sorted_reads_qc/stats/multiqc_report.html",
#         multiqc_DamID_reads="DamID_reads_bam_qc/stats/multiqc_report.html",
#         DamID_reads_bam_qc=expand("DamID_reads_bam_qc/stats/{sample}.stats.txt", sample=SAMPLES),

#         qc1="DamID_reads_bam_qc/fingerprint.pdf",
#         #qc2="DamID_reads_bam_qc/fragment_size.pdf",
#         qc3="DamID_reads_bam_qc/multiBamSummary.heatmap.pdf",
#         qc4="DamID_reads_bam_qc/multiBamSummary.pca.pdf",
#         qc5=expand("DamID_reads_bam_qc/{sample}.insert_size.pdf", sample=SAMPLES),
#         dag="Workflow_DAG.all.svg"







if DATA_TYPE == 'DamID' and config['MODE'] == 'SITE':
    rule macs2_DamID_sample_SITE:
        input:
            "results/clean_reads/{sample}.bam"
        output:
            "results/macs2_sample_level_peaks/{sample}_peaks.narrowPeak", 
            temp("results/macs2_sample_level_peaks/{sample}_treat_pileup.bdg"),
            temp("results/macs2_sample_level_peaks/{sample}_control_lambda.bdg")
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 8000  
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
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 8000
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
    rule macs2_sample_PE:
        input:
            "results/clean_reads/{sample}.bam"
        output:
            "results/macs2_sample_level_peaks/{sample}_peaks.narrowPeak", 
            temp("results/macs2_sample_level_peaks/{sample}_treat_pileup.bdg"),
            temp("results/macs2_sample_level_peaks/{sample}_control_lambda.bdg")
        params:
            odir='results/macs2_sample_level_peaks'
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 8000
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
             -n {wildcards.sample} --outdir {params.odir} -B &> {log}
            """

    rule macs2_contrast_PE:
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
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 8000
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
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1600
    threads:
        1
    log:
        "log/macs2_sample_level_peaks/{sample}_treat_pileup.bw.log"
    benchmark:
        "log/macs2_sample_level_peaks/{sample}_treat_pileup.bw.benchmark"
    priority:
        100
    conda:
        "../envs/bdg2bw.yaml"
    shell:
        """
        sort -k1,1 -k2,2n {input} > {output.sbdg} 2> {log}
        bedGraphToBigWig {output.sbdg} {SizeFile} {output.bw} &>> {log}
        """

# rule clean_peaks todo

rule peak2gtf:
    input:
        "results/macs2_sample_level_peaks/{sample}_peaks.narrowPeak"
    output:
        "results/macs2_sample_level_peaks/{sample}_peaks.gtf"
    threads:
        1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2000
    log:
        "log/macs2_sample_level_peaks_peak2gtf/{sample}.log"
    shell:
        """
        perl ../scripts/peak2gtf.pl {input} > {output} 
        """

rule peak_count:
    # todo broad/narrow
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
        maxFragmentLength=maxFragmentLength,
        minFragmentLength=minFragmentLength,
        MQ_MIN=MQ_MIN
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    conda: 
        '../envs/chiplike.yaml'
    shell:
        """
        featureCounts -a {input.gtf} -o {output} \
        -T {threads} -g gene_id -t peak -s 0 -p -B -C -d {params.minFragmentLength} -D {params.maxFragmentLength} \
        -Q {params.MQ_MIN} --minOverlap 1 --fracOverlap 0 \
        {input.bam} &> {log}
        """

    

rule bamCoverage_PE:
    # for ChIP
    input:
        "results/clean_reads/{sample}.bam"
    output:
        "results/clean_reads_bigWig/{sample}.cpm.bw"
    threads:
        8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    log:
        "log/clean_reads_bigWig/{sample}.bamCoverage.log"
    benchmark:
        "log/clean_reads_bigWig/{sample}.bamCoverage.benchmark"
    conda: 
        "../envs/deeptools.yaml"
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
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 6000
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
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
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
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
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
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
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
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000
    threads:
        1
    priority:
        100
    log:
        "log/macs2_contrast_level_peaks/{contrast}/log/{contrast_name}_control_lambda.bw.log"
    benchmark:
        "log/macs2_contrast_level_peaks/{contrast}/log/{contrast_name}_control_lambda.bw.benchmark"
    conda:
        "../envs/bdg2bw.yaml"
    shell:
        """
        sort -k1,1 -k2,2n {input} > {output.sbdg} 2> {log}
        bedGraphToBigWig {output.sbdg} \
            {SizeFile} {output.bw} &>> {log}
        """

rule contrast_treat_pileup_bw:
    input:
        "results/macs2_contrast_level_peaks/{contrast}/{contrast_name}_treat_pileup.bdg"
    output:
        bw="results/macs2_contrast_level_peaks/{contrast}/{contrast_name}_treat_pileup.bw",
        sbdg=temp("results/macs2_contrast_level_peaks/{contrast}/{contrast_name}_treat_pileup.s.bdg")
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000
    threads:
        1
    priority:
        100
    log:
        "log/macs2_contrast_level_peaks/{contrast}/{contrast_name}_treat_pileup.bw.log"
    benchmark:
        "log/macs2_contrast_level_peaks/{contrast}/{contrast_name}_treat_pileup.bw.benchmark"
    conda:
        "../envs/bdg2bw.yaml"
    shell:
        """
        sort -k1,1 -k2,2n {input} > {output.sbdg} 2>> {log}
        bedGraphToBigWig {output.sbdg} \
            {SizeFile} {output.bw} &>> {log}
        """



