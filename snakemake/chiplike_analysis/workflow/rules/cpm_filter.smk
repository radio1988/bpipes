
MIN_LFC=config['MIN_LFC']
MAX_FDR=config['MAX_FDR']


rule cpm_filter:
    input:
        COUNT="results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_count.txt",
        PEAK="results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.{narrowbroad}Peak",
    output:
        PEAK_UP="results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.real.{narrowbroad}Peak",
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
        Rscript workflow/scripts/cpm_filter.R {input.COUNT} {input.PEAK} {output.PEAK_UP} \
        config/meta.csv config/contrast.csv &> {log}
        """