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
    get_treat_pileup_bw_names_from_contrasts, \
    get_control_lambda_bw_names_from_contrasts


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
BW_BIN_SIZE=config['BW_BIN_SIZE']

#o=parse_meta_contrast(fmeta=pwd+"/config/meta.csv", fcontrast=pwd + "/config/contrast.csv") 
# print("parse_meta_contrast_obj:", vars(o))
# {'contrast2contrast_name': {'contrast1': 'G1_vs_ctrl', 'contrast2': 'G2_vs_ctrl', 'contrast3': 'G1_G2_vs_ctrl'}, 
# 'contrast2treatmentSamples': {'contrast1': ['2-1_S2', '2-2_S3', '2-3_S4'], 'contrast2': ['3-1_S5', '3-2_S6', '3-3_S7']}, 
# 'contrast2controlSamples': {'contrast1': ['1-2_S1'], 'contrast2': ['1-2_S1'], 'contrast3': ['1-2_S1']}}

ruleorder:  contrast_treat_pileup_bw > contrast_control_lambda_bw > sample_bdg2bw  
# Learn: To Avoid AmbiguousRuleException:
# Rules macs2_DamID_contrast_treat_pileup_bw and macs2_DamID_sample_treat_pileup_bw are ambiguous for the file narrow_peaks_sample_level/1-2_S1_treat_pileup.bw.
# Consider starting rule output with a unique prefix, constrain your wildcards, or use the ruleorder directive.
# Wildcards:
#     macs2_DamID_contrast_treat_pileup_bw: sample=1-2_S1
#     macs2_DamID_sample_treat_pileup_bw: sample=1-2_S1
# Expected input files:
#     macs2_DamID_contrast_treat_pileup_bw: narrow_peaks_sample_level/1-2_S1_treat_pileup.bdg
#     macs2_DamID_sample_treat_pileup_bw: narrow_peaks_sample_level/1-2_S1_treat_pileup.bdgExpected output files:
#     macs2_DamID_contrast_treat_pileup_bw: narrow_peaks_sample_level/1-2_S1_treat_pileup.bw
#     macs2_DamID_sample_treat_pileup_bw: narrow_peaks_sample_level/1-2_S1_treat_pileup.bw


    

rule sample_bdg2bw:
    input:
        "results/narrow_peaks_sample_level/{sample}_treat_pileup.bdg"
    output:
        bw="results/narrow_peaks_sample_level/{sample}_treat_pileup.bw",
        sbdg=temp("results/narrow_peaks_sample_level/{sample}_treat_pileup.s.bdg")
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1600
    threads:
        1
    log:
        "log/narrow_peaks_sample_level/{sample}_treat_pileup.bw.log"
    benchmark:
        "log/narrow_peaks_sample_level/{sample}_treat_pileup.bw.benchmark"
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

    
rule bamCoverage:
    # for ChIP
    input:
        "results/clean_reads/{sample}.bam"
    output:
        "results/clean_reads_bigWig/{sample}.cpm.bw"
    params:
        bin_size=BW_BIN_SIZE,
        pse=("--minFragmentLength {config[minFragmentLength]} --maxFragmentLength {config[maxFragmentLength]}  -e 150"
                if MODE=='PE'
                else "" )
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
        --binSize {params.bin_size} \
        {params.pse} \
        &> {log}
        """

### narrow_peaks_contrast_level ### 


rule contrast_control_lambda_bw:
    input:
        "results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_control_lambda.bdg"
    output:
        bw="results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_control_lambda.bw",
        sbdg=temp("results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_control_lambda.s.bdg")
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000
    threads:
        1
    priority:
        100
    log:
        "log/narrow_peaks_contrast_level/{contrast}/log/{contrast_name}_control_lambda.bw.log"
    benchmark:
        "log/narrow_peaks_contrast_level/{contrast}/log/{contrast_name}_control_lambda.bw.benchmark"
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
        "results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_treat_pileup.bdg"
    output:
        bw="results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_treat_pileup.bw",
        sbdg=temp("results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_treat_pileup.s.bdg")
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000
    threads:
        1
    priority:
        100
    log:
        "log/narrow_peaks_contrast_level/{contrast}/{contrast_name}_treat_pileup.bw.log"
    benchmark:
        "log/narrow_peaks_contrast_level/{contrast}/{contrast_name}_treat_pileup.bw.benchmark"
    conda:
        "../envs/bdg2bw.yaml"
    shell:
        """
        sort -k1,1 -k2,2n {input} > {output.sbdg} 2>> {log}
        bedGraphToBigWig {output.sbdg} \
            {SizeFile} {output.bw} &>> {log}
        """



