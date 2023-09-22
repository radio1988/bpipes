BLACKLIST=config['BLACKLIST']
MODE=config['MODE']
GSIZE=config['GSIZE']
SizeFile=config['SizeFile']
MAX_FDR=config['MAX_FDR']



if DATA_TYPE == 'DamID' and 'MODE' == 'SITE':
    rule macs2_DamID_sample_SITE:
        input:
            "results/clean_reads/{sample}.bam"
        output:
            "results/narrow_peaks_sample_level/{sample}_peaks.narrowPeak", 
            temp("results/narrow_peaks_sample_level/{sample}_treat_pileup.bdg"),
            temp("results/narrow_peaks_sample_level/{sample}_control_lambda.bdg")
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 8000  
        threads:
            1
        conda:
            "../envs/macs2.yaml"
        log:
            "log/narrow_peaks_sample_level/{sample}_peaks.narrowPeak.log"
        benchmark:
            "log/narrow_peaks_sample_level/{sample}_peaks.narrowPeak.benchmark"
        shell:
            """
             macs2 callpeak -t {input} \
             -f BAM --nomodel --shift -60 --extsize 100 -g {GSIZE} -q {MAX_FDR} --keep-dup all \
             -n {wildcards.sample} --outdir results/narrow_peaks_sample_level -B &> {log}
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
            "results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_peaks.narrowPeak", # e.g. "narrow_peaks_contrast_level/contrast1/G1_vs_ctrl_peaks.narrowPeak"
            "results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_summits.bed", 
            temp("results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_treat_pileup.bdg"),
            temp("results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_control_lambda.bdg")
        params:
            contrast_name=lambda wildcards: get_contrast_name_from_contrast(contrast=wildcards.contrast, o=o),
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 8000
        threads:
            1        
        conda:
            "../envs/macs2.yaml"
        log:
            "log/narrow_peaks_contrast_level/{contrast}/{contrast_name}.macs2_DamID.log"
        benchmark:
            "log/narrow_peaks_contrast_level/{contrast}/{contrast_name}.macs2_DamID.benchmark"
        shell:
            """
            macs2 callpeak -t {input.treatment} -c {input.control} \
            -f BAM --nomodel --shift -60 --extsize 100 -g {GSIZE} -q {MAX_FDR} --keep-dup all \
            -n {params.contrast_name} --outdir results/narrow_peaks_contrast_level/{wildcards.contrast} -B &> {log}
            """
elif (DATA_TYPE == 'DamID' or DATA_TYPE=='ChIP') and MODE in ['PE', 'SE']:
    rule macs2_narrow_sample:
        input:
            "results/clean_reads/{sample}.bam"
        output:
            "results/narrow_peaks_sample_level/{sample}/{sample}_peaks.narrowPeak", 
            temp("results/narrow_peaks_sample_level/{sample}/{sample}_treat_pileup.bdg"),
            temp("results/narrow_peaks_sample_level/{sample}/{sample}_control_lambda.bdg")
        params:
            odir='results/narrow_peaks_sample_level/{sample}',
            pse='BAMPE'
                    if MODE=='PE'
                    else 'BAM'
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 8000
        threads:
            1
        conda:
            "../envs/macs2.yaml"
        log:
            "log/narrow_peaks_sample_level/{sample}_peaks.narrowPeak.log"
        benchmark:
            "log/narrow_peaks_sample_level/{sample}_peaks.narrowPeak.benchmark"
        shell:
            """
             macs2 callpeak -t {input} \
             -f {params.pse} -g {GSIZE} -q {MAX_FDR} --keep-dup all \
             -n {wildcards.sample} --outdir {params.odir} -B &> {log}
            """

    rule macs2_broad_sample:
        input:
            "results/clean_reads/{sample}.bam"
        output:
            "results/broad_peaks_sample_level/{sample}/{sample}_peaks.broadPeak", 
        params:
            odir='results/broad_peaks_sample_level/{sample}',
            pse='BAMPE'
                    if MODE=='PE'
                    else 'BAM'
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 8000
        threads:
            1
        conda:
            "../envs/macs2.yaml"
        log:
            "log/broad_peaks_sample_level/{sample}_peaks.broadPeak.log"
        benchmark:
            "log/broad_peaks_sample_level/{sample}_peaks.broadPeak.benchmark"
        shell:
            """
             macs2 callpeak -t {input} \
             -f {params.pse} -g {GSIZE} -q {MAX_FDR} --keep-dup all \
             -n {wildcards.sample} --outdir {params.odir} \
             --broad --broad-cutoff 0.1 &> {log}
            """

    rule macs2_narrow_contrast:
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
            "results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_peaks.narrowPeak",
            temp("results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_summits.bed"), # for narrowPeak
            temp("results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_treat_pileup.bdg"),
            temp("results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_control_lambda.bdg")
        params:
            contrast_name=lambda wildcards: get_contrast_name_from_contrast(contrast=wildcards.contrast, o=o),
            pse='BAMPE'
                    if MODE=='PE'
                    else 'BAM'
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 8000
        threads:
            1        
        conda:
            "../envs/macs2.yaml"
        log:
            "log/narrow_peaks_contrast_level/{contrast}/{contrast_name}.log"
        benchmark:
            "log/narrow_peaks_contrast_level/{contrast}/{contrast_name}.benchmark"
        shell:
            """
            macs2 callpeak -t {input.treatment} -c {input.control} \
            -f {params.pse} -g {GSIZE} -q {MAX_FDR} --keep-dup all \
            -n {params.contrast_name} --outdir results/narrow_peaks_contrast_level/{wildcards.contrast} -B &> {log}
            """

    rule macs2_broad_contrast:
        input:
            treatment=lambda wildcards: get_treatment_bams_from_contrast(contrast=wildcards.contrast, o=o),
            control=lambda wildcards: get_control_bams_from_contrast(contrast=wildcards.contrast, o=o),
        output:
            "results/broad_peaks_contrast_level/{contrast}/{contrast_name}_peaks.broadPeak",
            # e.g. "broad_peaks_contrast_level/contrast1/G1_vs_ctrl_peaks.broadPeak"
            temp("results/broad_peaks_contrast_level/{contrast}/{contrast_name}_peaks.gappedPeak"), 
            temp("results/broad_peaks_contrast_level/{contrast}/{contrast_name}_peaks.xls"), 
        params:
            contrast_name=lambda wildcards: get_contrast_name_from_contrast(contrast=wildcards.contrast, o=o),
            pse='BAMPE'
                    if MODE=='PE'
                    else 'BAM'
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 8000
        threads:
            1        
        conda:
            "../envs/macs2.yaml"
        log:
            "log/broad_peaks_contrast_level/{contrast}/{contrast_name}.log"
        benchmark:
            "log/broad_peaks_contrast_level/{contrast}/{contrast_name}.benchmark"
        shell:
            """
            macs2 callpeak -t {input.treatment} -c {input.control} \
            -f {params.pse} -g {GSIZE} -q {MAX_FDR} --keep-dup all \
            -n {params.contrast_name} --outdir results/broad_peaks_contrast_level/{wildcards.contrast} \
            --broad --broad-cutoff 0.1 &> {log}
            """
else:
    sys.exit("MODE Error")



rule blacklist_filter:
    input:
        peak="results/{narrowbroad}_peaks_{cs}_level/{name}/{name2}_peaks.{narrowbroad}Peak",
        blacklist=BLACKLIST
    output:
        "results/{narrowbroad}_peaks_{cs}_level/{name}/{name2}_clean.{narrowbroad}Peak"
    log:
        "log/results/{narrowbroad}_peaks_{cs}_level/{name}/{name2}_clean.{narrowbroad}Peak.log"
    benchmark:
        "log/results/{narrowbroad}_peaks_{cs}_level/{name}/{name2}_clean.{narrowbroad}Peak.benchmark"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    threads:
        1
    conda:
        "../envs/bedtools.yaml"   
    shell:
        """
        bedtools subtract -v -a {input.peak} -b {input.blacklist} > {output}
        """


### Tracks bdg2bw

ruleorder:  contrast_treat_bdg2bw > contrast_control_bdg2bw > sample_bdg2bw
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
#     macs2_DamID_sample_tre

rule contrast_control_bdg2bw:
    input:
        "results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_control_lambda.bdg"
    output:
        bw="results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_control_lambda.bw",
        sbdg=temp("results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_control_lambda.s.bdg")
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 32000
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
        LC_COLLATE=C sort -k1,1 -k2,2n {input} > {output.sbdg} 2> {log}
        bedGraphToBigWig {output.sbdg} \
            {SizeFile} {output.bw} &>> {log}
        """

rule contrast_treat_bdg2bw:
    input:
        "results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_treat_pileup.bdg"
    output:
        bw="results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_treat_pileup.bw",
        sbdg=temp("results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_treat_pileup.s.bdg")
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 32000
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
        LC_COLLATE=C sort -k1,1 -k2,2n {input} > {output.sbdg} 2>> {log}
        bedGraphToBigWig {output.sbdg} \
            {SizeFile} {output.bw} &>> {log}
        """

rule sample_bdg2bw:
    input:
        "results/narrow_peaks_sample_level/{sample}_treat_pileup.bdg"
    output:
        bw="results/narrow_peaks_sample_level/{sample}_treat_pileup.bw",
        sbdg=temp("results/narrow_peaks_sample_level/{sample}_treat_pileup.s.bdg")
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 32000
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
        LC_COLLATE=C sort -k1,1 -k2,2n {input} > {output.sbdg} 2> {log}
        bedGraphToBigWig {output.sbdg} {SizeFile} {output.bw} &>> {log}
        """

