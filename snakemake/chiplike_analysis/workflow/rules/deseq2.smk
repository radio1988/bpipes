
MIN_LFC=config['MIN_LFC']
MAX_FDR=config['MAX_FDR']


rule deseq2:
    input:
        count="results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_count.txt",
        anno="results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.{narrowbroad}Peak.full_anno.xlsx",
        config='config/config.yaml' # update if this update
    output:
        full="results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast}.{contrast_name}.DESeq2.xlsx",
        sig=
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000
    threads:
        1
    log:
        "log/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast}.{contrast_name}.DESeq2.log"
    benchmark:
        "log/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast}.{contrast_name}.DESeq2.benchmark"
    conda:
        "../envs/deseq2.yaml"
    shell:
        """
        Rscript workflow/scripts/deseq2.R {input.count} {input.anno} \
        config/meta.csv config/contrast.csv {MAX_FDR} {MIN_LFC} &> {log}
        """


rule deseq2_peak_filter:
    input:
        peak="results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.{narrowbroad}Peak",
        sig="results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast}.{contrast_name}.DESeq2.xlsx",
        config='config/config.yaml'
    output:
        "results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.sig.{narrowbroad}Peak"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000
    threads:
        1
    log:
        "log/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.sig.{narrowbroad}Peak.log"
    conda:
        "../envs/deseq2.yaml"
    shell:
        """
        Rscript workflow/scripts/deseq2.filter_peaks.R {input.peak} {input.sig} {MAX_FDR} {MIN_LFC} &> {log}
        """

