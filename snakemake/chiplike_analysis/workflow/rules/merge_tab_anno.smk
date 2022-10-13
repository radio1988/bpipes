GTF=config['GTF']
BIDING_LEFT=config['BIDING_LEFT']
BIDING_RIGHT=config['BIDING_RIGHT']
CHIPPEAKANNO_MODE=config['CHIPPEAKANNO_MODE']


rule merge_tab_anno:
    input:
        ANNO="results/{narrowbroad}_peaks_contrast_level/" \
             + "{contrast}/{contrast_name}_clean.{narrowbroad}Peak.full_anno.xlsx",
        TAB="results/{narrowbroad}_peaks_contrast_level/"\
            + "{contrast}/{contrast_name}_clean.t_vs_c.xlsx"
    output:
        "results/{narrowbroad}_peaks_contrast_level/" \
             + "{contrast}/{contrast_name}_clean.{narrowbroad}Peak.final_anno.xlsx"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000
    threads:
        1
    log:
        "log/{narrowbroad}_peaks_contrast_level/" \
             + "{contrast}/{contrast_name}_clean.{narrowbroad}Peak.final_anno.log"
    benchmark:
        "log/{narrowbroad}_peaks_contrast_level/" \
             + "{contrast}/{contrast_name}_clean.{narrowbroad}Peak.final_anno.benchmark"
    conda:
        "../envs/deseq2.yaml"
    shell:
        """
        Rscript workflow/scripts/merge_table.R {input.ANNO} {input.TAB} {output} &> {log}
        """