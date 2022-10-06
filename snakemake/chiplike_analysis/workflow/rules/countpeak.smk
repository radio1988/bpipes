SAMPLES=config['SAMPLES']

rule peak2gtf_sample_level:
    input:
        "results/{narrowbroad}_peaks_sample_level/{name1}/{name2}_clean.{narrowbroad}Peak"
    output:
        "results/{narrowbroad}_peaks_sample_level/{name1}/{name2}_clean.gtf"
    threads:
        1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2000
    conda: 
        '../envs/chiplike.yaml'
    log:
        "log/{narrowbroad}_peaks_sample_level_peak2gtf/{name1}/{name2}.log"
    shell:
        """
        perl workflow/scripts/peak2gtf.pl {input} > {output} 2> {log}
        """

rule peak_count_sample_level_pe:
    # todo broad/narrow
    input:
        bam=expand("results/clean_reads/{sample}.bam", sample=SAMPLES),
        gtf="results/{narrowbroad}_peaks_sample_level/{name1}/{name2}_clean.gtf"
    output:
        "results/{narrowbroad}_peaks_sample_level/{name1}/{name2}.count.txt"
    log:
        "log/{narrowbroad}_peaks_sample_level/{name1}/{name2}.count.log"
    benchmark:
        "log/{narrowbroad}_peaks_sample_level/{name1}/{name2}.count.benchmark"
    threads:
        4
    params:
        maxFragmentLength=maxFragmentLength,
        minFragmentLength=minFragmentLength,
        MQ_MIN=MQ_MIN
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    conda: 
        '../envs/subread.yaml'
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
        "results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.{narrowbroad}Peak"
    output:
        "results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.gtf"
    threads:
        1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2000
    conda: 
        '../envs/chiplike.yaml'
    log:
        "log/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.gtf.log"
    shell:
        """
        perl workflow/scripts/peak2gtf.pl {input} > {output} 2> {log}
        """

rule peak_count_contrast_level_pe:
    # todo broad/narrow
    input:
        bam=expand("results/clean_reads/{sample}.bam", sample=SAMPLES),
        gtf="results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.gtf"
    output:
        "results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_count.txt"
    log:
        "log/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_count.txt.log"
    benchmark:
        "log/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_count.txt.log"
    threads:
        4
    params:
        maxFragmentLength=maxFragmentLength,
        minFragmentLength=minFragmentLength,
        MQ_MIN=MQ_MIN
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    conda: 
        '../envs/subread.yaml'
    shell:
        """
        featureCounts -a {input.gtf} -o {output} \
        -T {threads} -g gene_id -t peak -s 0 -p -B -C -d {params.minFragmentLength} -D {params.maxFragmentLength} \
        -Q {params.MQ_MIN} --minOverlap 1 --fracOverlap 0 \
        {input.bam} &> {log}
        """
