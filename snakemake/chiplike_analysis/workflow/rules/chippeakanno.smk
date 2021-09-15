GTF=config['GTF']
BIDING_LEFT=config['BIDING_LEFT']
BIDING_RIGHT=config['BIDING_RIGHT']
CHIPPEAKANNO_MODE=config['CHIPPEAKANNO_MODE']


rule chippeakanno:
    input:
        peak="results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.real.{narrowbroad}Peak",
        gtf=GTF,
        config='config/config.yaml'
    output:
        "results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.{narrowbroad}Peak.full_anno.xlsx"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000
    threads:
        1
    log:
        "log/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.{narrowbroad}Peak.full_anno.xlsx.log"
    benchmark:
        "log/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.{narrowbroad}Peak.full_anno.xlsx.benchmark"
    conda:
        "../envs/chippeakanno.yaml"
    shell:
        """
        Rscript workflow/scripts/chippeakanno.R {input.peak} {input.gtf} \
        {CHIPPEAKANNO_MODE} {BIDING_LEFT} {BIDING_RIGHT} &> {log}
        """