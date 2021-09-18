MIN_LFC=config['MIN_LFC']
MAX_FDR=config['MAX_FDR']


rule t_vs_c_peak_filter:
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



rule t_vs_c_summit_filter:
    '''
    treatment (Pulldown) should be higher than control (input/IgG) in normazlied count
    output called "real" peaks
    '''
    input:
        PEAK_UP="results/narrow_peaks_contrast_level/"\
            +"{contrast}/{contrast_name}_clean.real.narrowPeak",
        SUMMIT="results/narrow_peaks_contrast_level/"\
             +"{contrast}/{contrast_name}_summits.bed",
    output:
        "results/narrow_peaks_contrast_level/"\
            +"{contrast}/{contrast_name}_summits.real.bed"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000
    threads:
        1
    log:
        "log/narrow_peaks_contrast_level/"\
            +"{contrast}/{contrast_name}_summits.real.bed.log"
    benchmark:
        "log/narrow_peaks_contrast_level/"\
            +"{contrast}/{contrast_name}_summits.real.bed.benchmark"
    conda:
        "../envs/deseq2.yaml"
    shell:
        """
        Rscript workflow/scripts/summit_filter.R {input.SUMMIT} {input.PEAK_UP} \
        {output} &> {log}
        """