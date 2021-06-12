SAMPLES=config['SAMPLES']

rule peak2gtf_sample_level:
    input:
        "results/narrow_peaks_sample_level/{sample}_clean.narrowPeak"
    output:
        "results/narrow_peaks_sample_level/{sample}_clean.gtf"
    threads:
        1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2000
    conda: 
        '../envs/chiplike.yaml'
    log:
        "log/narrow_peaks_sample_level_peak2gtf/{sample}.log"
    shell:
        """
        perl workflow/scripts/peak2gtf.pl {input} > {output} 2> {log}
        """

rule peak_count_sample_level:
    # todo broad/narrow
    input:
        bam=expand("results/clean_reads/{sample}.bam", sample=SAMPLES),
        gtf="results/narrow_peaks_sample_level/{sample}_clean.gtf"
    output:
        "results/narrow_peaks_sample_level/{sample}.count.txt"
    log:
        "log/narrow_peaks_sample_level/{sample}.count.log"
    benchmark:
        "log/narrow_peaks_sample_level/{sample}.count.benchmark"
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

rule peak2gtf_contrast_level:
# todo: handle zero peak gtf files
    input:
        "results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_clean.narrowPeak"
    output:
        "results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_clean.gtf"
    threads:
        1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2000
    conda: 
        '../envs/chiplike.yaml'
    log:
        "log/narrow_peaks_contrast_level/{contrast}/{contrast_name}_clean.gtf.log"
    shell:
        """
        perl workflow/scripts/peak2gtf.pl {input} > {output} 2> {log}
        """

rule peak_count_contrast_level:
    # todo broad/narrow
    input:
        bam=expand("results/clean_reads/{sample}.bam", sample=SAMPLES),
        gtf="results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_clean.gtf"
    output:
        "results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_count.txt"
    log:
        "log/narrow_peaks_contrast_level/{contrast}/{contrast_name}_count.txt.log"
    benchmark:
        "log/narrow_peaks_contrast_level/{contrast}/{contrast_name}_count.txt.log"
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