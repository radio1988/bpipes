SAMPLES=config['SAMPLES']
GTF=config['GTF']
signalHeatmapWidth=config['signalHeatmapWidth']

def clean_or_sorted_bams_input(wildcards):
   return ["results/"+wildcards.cs_folder+"/"+sample+".bam" \
            for sample in SAMPLES]


## TSS based
rule signalHeatmapMatrix_TSS:
    input:
        BWS=expand("results/clean_reads_bigWig/{sample}.cpm.bw", sample=SAMPLES),
        gtf=GTF
    output:
        "results/clean_reads_qc/signalHeatmap/clean.all_sample.mat.gz"
    params:
        a=signalHeatmapWidth, # downstream
        b=signalHeatmapWidth  # upstream
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000
    threads:
        4
    log:
        "log/clean_reads_qc/signalHeatmap/clean.all_sample.mat.gz.log"
    benchmark:
        "log/clean_reads_qc/signalHeatmap/clean.all_sample.mat.gz.benchmark"
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        computeMatrix reference-point -S {input.BWS} \
            -R {input.gtf} \
            --referencePoint TSS -a {params.a} -b {params.b} \
            -p {threads} \
            -out {output}  &> {log}
        """

   
rule signalHeatmap_TSS:
    """
    colormap: https://matplotlib.org/2.0.2/users/colormaps.html
    jet (rainbow), viridis, seismic (blue-red)
    """
    input:
        "results/clean_reads_qc/signalHeatmap/clean.all_sample.mat.gz"
    output:
        "results/clean_reads_qc/signalHeatmap/clean.all_sample.TSS.signalHeatmap.pdf"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000
    threads:
        4
    log:
        "log/clean_reads_qc/signalHeatmap/clean.all_sample.TSS.signalHeatmap.pdf.log"
    benchmark:
        "log/clean_reads_qc/signalHeatmap/clean.all_sample.TSS.signalHeatmap.pdf.benchmark"
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        plotHeatmap -m {input} \
            --zMin 0 --zMax auto --colorMap viridis \
            --refPointLabel 'TSS' \
            -T 'Signal Around TSS' \
            -y 'signal' -x '' -z 'transcripts' \
            --verbose  -out {output} &> {log}
        """

## Peak based
rule signalHeatmapMatrix_ContrastPeak:
    input:
        BWS=expand("results/clean_reads_bigWig/{sample}.cpm.bw", sample=SAMPLES),
        peak="results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.{narrowbroad}Peak",
    output:
        "results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.{narrowbroad}Peak.mat.gz"
    params:
        a=signalHeatmapWidth, # downstream
        b=signalHeatmapWidth  # upstream
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 64000
    threads:
        4
    log:
        "log/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.{narrowbroad}Peak.mat.gz.log"
    benchmark:
        "log/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.{narrowbroad}Peak.mat.gz.benchmark"
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        computeMatrix reference-point \
            -S {input.BWS} \
            -R {input.peak} \
            --referencePoint center -a {params.a} -b {params.b} \
            -p {threads} \
            -out {output}  &> {log}
        """

rule signalHeatmap_ContrastPeak:
    input:
        "results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.{narrowbroad}Peak.mat.gz"
    output:
        "results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.{narrowbroad}Peak.pdf"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 64000
    threads:
        1
    log:
        "log/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.{narrowbroad}Peak.pdf.log"
    benchmark:
        "log/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.{narrowbroad}Peak.pdf.benchmark"
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        plotHeatmap -m {input} \
            --zMin 0 --zMax auto --colorMap viridis \
            --refPointLabel 'peak_center' \
            -T 'Signal Around Peak Center' \
            -y 'signal' -x '' -z 'peaks' \
            --verbose -out {output} &> {log}
        """

