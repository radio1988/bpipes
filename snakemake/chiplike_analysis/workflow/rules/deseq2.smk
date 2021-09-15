
MIN_LFC=config['MIN_LFC']
MAX_FDR=config['MAX_FDR']


rule deseq2:
    input:
        COUNT="results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_count.txt",
        ANNO="results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.{narrowbroad}Peak.full_anno.xlsx",
        config='config/config.yaml' # update if this update
    output:
        "results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast}.{contrast_name}.DESeq2.xlsx"
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
        Rscript workflow/scripts/deseq2.R {input.COUNT} {input.ANNO} \
        config/meta.csv config/contrast.csv {MAX_FDR} {MIN_LFC} &> {log}
        """