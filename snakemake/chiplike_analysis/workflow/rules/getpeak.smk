BLACKLIST=config['BLACKLIST']




if DATA_TYPE == 'DamID' and config['MODE'] == 'SITE':
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
            4
        conda:
            "../envs/macs2.yaml"
        log:
            "log/narrow_peaks_sample_level/{sample}_peaks.narrowPeak.log"
        benchmark:
            "log/narrow_peaks_sample_level/{sample}_peaks.narrowPeak.benchmark"
        shell:
            """
             macs2 callpeak -t {input} \
             -f BAM --nomodel --shift -60 --extsize 100 -g {GSIZE} -q 0.05 --keep-dup all \
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
            4        
        conda:
            "../envs/macs2.yaml"
        log:
            "log/narrow_peaks_contrast_level/{contrast}/{contrast_name}.macs2_DamID.log"
        benchmark:
            "log/narrow_peaks_contrast_level/{contrast}/{contrast_name}.macs2_DamID.benchmark"
        shell:
            """
            macs2 callpeak -t {input.treatment} -c {input.control} \
            -f BAM --nomodel --shift -60 --extsize 100 -g {GSIZE} -q 0.05 --keep-dup all \
            -n {params.contrast_name} --outdir results/narrow_peaks_contrast_level/{wildcards.contrast} -B &> {log}
            """
elif (DATA_TYPE == 'DamID' or DATA_TYPE=='ChIP') and config['MODE'] == 'PE':
    rule macs2_narrow_sample_pe:
        input:
            "results/clean_reads/{sample}.bam"
        output:
            "results/narrow_peaks_sample_level/{sample}/{sample}_peaks.narrowPeak", 
            temp("results/narrow_peaks_sample_level/{sample}/{sample}_treat_pileup.bdg"),
            temp("results/narrow_peaks_sample_level/{sample}/{sample}_control_lambda.bdg")
        params:
            odir='results/narrow_peaks_sample_level/{sample}'
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 8000
        threads:
            4
        conda:
            "../envs/macs2.yaml"
        log:
            "log/narrow_peaks_sample_level/{sample}_peaks.narrowPeak.log"
        benchmark:
            "log/narrow_peaks_sample_level/{sample}_peaks.narrowPeak.benchmark"
        shell:
            """
             macs2 callpeak -t {input} \
             -f BAMPE -g {GSIZE} -q 0.05 --keep-dup all \
             -n {wildcards.sample} --outdir {params.odir} -B &> {log}
            """

    rule macs2_broad_sample_pe:
        input:
            "results/clean_reads/{sample}.bam"
        output:
            "results/broad_peaks_sample_level/{sample}/{sample}_peaks.broadPeak", 
        params:
            odir='results/broad_peaks_sample_level/{sample}'
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 8000
        threads:
            4
        conda:
            "../envs/macs2.yaml"
        log:
            "log/broad_peaks_sample_level/{sample}_peaks.broadPeak.log"
        benchmark:
            "log/broad_peaks_sample_level/{sample}_peaks.broadPeak.benchmark"
        shell:
            """
             macs2 callpeak -t {input} \
             -f BAMPE -g {GSIZE} -q 0.05 --keep-dup all \
             -n {wildcards.sample} --outdir {params.odir} \
             --broad --broad-cutoff 0.1 &> {log}
            """

    rule macs2_narrow_contrast_pe:
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
            "results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_summits.bed", # for narrowPeak
            temp("results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_treat_pileup.bdg"),
            temp("results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_control_lambda.bdg")
        params:
            contrast_name=lambda wildcards: get_contrast_name_from_contrast(contrast=wildcards.contrast, o=o),
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 8000
        threads:
            4        
        conda:
            "../envs/macs2.yaml"
        log:
            "log/narrow_peaks_contrast_level/{contrast}/{contrast_name}.log"
        benchmark:
            "log/narrow_peaks_contrast_level/{contrast}/{contrast_name}.benchmark"
        shell:
            """
            macs2 callpeak -t {input.treatment} -c {input.control} \
            -f BAMPE -g {GSIZE} -q 0.05 --keep-dup all \
            -n {params.contrast_name} --outdir results/narrow_peaks_contrast_level/{wildcards.contrast} -B &> {log}
            """

    rule macs2_broad_contrast_pe:
        input:
            treatment=lambda wildcards: get_treatment_bams_from_contrast(contrast=wildcards.contrast, o=o),
            control=lambda wildcards: get_control_bams_from_contrast(contrast=wildcards.contrast, o=o),
        output:
            "results/broad_peaks_contrast_level/{contrast}/{contrast_name}_peaks.broadPeak", 
            # e.g. "broad_peaks_contrast_level/contrast1/G1_vs_ctrl_peaks.broadPeak"
            "results/broad_peaks_contrast_level/{contrast}/{contrast_name}_peaks.gappedPeak", 
            "results/broad_peaks_contrast_level/{contrast}/{contrast_name}_peaks.xls", 
        params:
            contrast_name=lambda wildcards: get_contrast_name_from_contrast(contrast=wildcards.contrast, o=o),
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 8000
        threads:
            4        
        conda:
            "../envs/macs2.yaml"
        log:
            "log/broad_peaks_contrast_level/{contrast}/{contrast_name}.log"
        benchmark:
            "log/broad_peaks_contrast_level/{contrast}/{contrast_name}.benchmark"
        shell:
            """
            macs2 callpeak -t {input.treatment} -c {input.control} \
            -f BAMPE -g {GSIZE} -q 0.05 --keep-dup all \
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
        bedtools intersect -v -a {input.peak} -b {input.blacklist} > {output}
        """


# rule blacklist_filter_contrast_level:
#     input:
#         peak="results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_peaks.{narrowbroad}Peak",
#         blacklist=BLACKLIST
#     output:
#         "results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.{narrowbroad}Peak"
#         # e.g. xxx_clean.narrowPeak, xxx_clean.broadPeak
#     log:
#         "log/results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.{narrowbroad}Peak.log"
#     benchmark:
#         "log/results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.{narrowbroad}Peak.benchmark"
#     resources:
#         mem_mb=lambda wildcards, attempt: attempt * 8000
#     threads:
#         1
#     conda:
#         "../envs/bedtools.yaml"   
#     shell:
#         """
#         bedtools intersect -v -a {input.peak} -b {input.blacklist} > {output}
#         """


# rule blacklist_filter_sample_level:
#     input:
#         peak="results/narrow_peaks_sample_level/{sample}_peaks.narrowPeak",
#         blacklist=BLACKLIST
#     output:
#         "results/narrow_peaks_sample_level/{sample}_clean.narrowPeak",
#     log:
#         "log/narrow_peaks_sample_level/{sample}_clean.narrowPeak.log",
#     benchmark:
#         "log/results/narrow_peaks_sample_level/{sample}_clean.narrowPeak.benchmark",
#     resources:
#         mem_mb=lambda wildcards, attempt: attempt * 8000
#     threads:
#         1
#     conda:
#         "../envs/bedtools.yaml"   
#     shell:
#         """
#         bedtools intersect -v -a {input.peak} -b {input.blacklist} > {output}
#         """











