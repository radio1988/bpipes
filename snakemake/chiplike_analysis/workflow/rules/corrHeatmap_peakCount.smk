rule corrHeatmap_peakCount:
    input:
        "results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_count.txt"
    output:
        "results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_count.heatmap.pdf"
    threads:
        1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2000
    conda: 
        '../envs/chiplike.yaml'
    log:
        "log/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_count.heatmap.pdf.log"
    benchmark:
        "log/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_count.heatmap.pdf.benchmark"
    shell:
        """
        Rscript workflow/scripts/corrHeatmap.peakCount.R {input} {output} &> {log}
        """


# todo: rule DESeq2
