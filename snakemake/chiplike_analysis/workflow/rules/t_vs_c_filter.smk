MIN_LFC=config['MIN_LFC']
MAX_FDR=config['MAX_FDR']


rule t_vs_c_filter:
    '''
    treatment (Pulldown) should be higher than control (input/IgG) in normazlied count
    output called "real" peaks
    '''
    input:
        COUNT="results/{narrowbroad}_peaks_contrast_level/"\
             +"{contrast}/{contrast_name}_count.txt",
        PEAK="results/{narrowbroad}_peaks_contrast_level/"\
             +"{contrast}/{contrast_name}_clean.{narrowbroad}Peak",
    output:
        PEAK_UP="results/{narrowbroad}_peaks_contrast_level/"\
                +"{contrast}/{contrast_name}_clean.real.{narrowbroad}Peak",
        TAB="results/{narrowbroad}_peaks_contrast_level/"\
                +"{contrast}/{contrast_name}_clean.t_vs_c.xlsx",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000
    threads:
        1
    log:
        "log/{narrowbroad}_peaks_contrast_level/"\
        +"{contrast}/{contrast_name}.t_vs_c_filter.log"
    benchmark:
        "log/{narrowbroad}_peaks_contrast_level/"\
        +"{contrast}/{contrast_name}.t_vs_c_filter.benchmark"
    conda:
        "../envs/deseq2.yaml"
    shell:
        """
        Rscript workflow/scripts/t_vs_c_filter.R {input.COUNT} {input.PEAK} {output.PEAK_UP} \
        config/meta.csv config/contrast.csv {output.TAB} &> {log}
        """