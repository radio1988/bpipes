GTF=config['GTF']
BIDING_LEFT=config['BIDING_LEFT']
BIDING_RIGHT=config['BIDING_RIGHT']
CHIPPEAKANNO_MODE=config['CHIPPEAKANNO_MODE']

def clean_or_sorted_bams_input(wildcards):
   return ["results/"+wildcards.cs_folder+"/"+sample+".bam" \
            for sample in SAMPLES]


rule chippeakanno:
    input:
        peak="results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.{narrowbroad}Peak",
        gtf=GTF
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
        Rscript workflow/scripts/chippeakanno.R {input.peak} {input.gtf} CHIPPEAKANNO_MODE BIDING_LEFT BIDING_RIGHT &> {log}
        """