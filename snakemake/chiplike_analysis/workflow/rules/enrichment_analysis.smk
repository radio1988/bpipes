ORG_EG_DB = config['ORG_EG_DB']


rule enrichment_analysis:
    '''
    treatment (Pulldown) should be higher than control (input/IgG) in normazlied count
    output called "real" peaks
    '''
    input:
        "results/{narrowbroad}_peaks_contrast_level/" \
             + "{contrast}/{contrast_name}_clean.{narrowbroad}Peak.final_anno.xlsx"
    params:
        ORG_EG_DB
    output:
        touch("results/{narrowbroad}_peaks_contrast_level/" \
             + "{contrast}/{contrast_name}_clean.real.{narrowbroad}Peak.enrichment.finished")
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000
    threads:
        1
    log:
        "log/{narrowbroad}_peaks_contrast_level/" \
             + "{contrast}/{contrast_name}_clean.real.{narrowbroad}Peak.enrichment.log"
    benchmark:
        "log/{narrowbroad}_peaks_contrast_level/" \
             + "{contrast}/{contrast_name}_clean.real.{narrowbroad}Peak.enrichment.benchmark"
    conda:
        "../envs/chippeakanno.yaml"
    shell:
        """
        Rscript workflow/scripts/enrichment_analysis.R {input} {params} &> {log}
        """
