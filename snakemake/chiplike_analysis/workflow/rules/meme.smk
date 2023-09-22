MEME_DB=config['MEME_DB']

ruleorder: get_peak_fasta_contrast_level > meme_peak_contrast_level > get_summit_neighbour > split_fa_by_chr > meme_neibour > meme_neibour_chr_split

## PEAK REGION ##

rule get_peak_fasta_contrast_level:
    input:
        peak="results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.{narrowbroad}Peak",
        genome=GENOME
    output:
        "results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.{narrowbroad}Peak.fa"
    log:
        "log/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.{narrowbroad}Peak.fa.log"
    threads:
        1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 6000
    conda:
        "../envs/biopython.yaml"
    shell:
        """
        python workflow/scripts/get_peak_fasta.py {input.genome} {input.peak} {output} &> {log}
        """

rule meme_peak_contrast_level:
    input: 
        fasta="results/{narrowbroad}_peaks_contrast_level/{contrast}/{contrast_name}_clean.{narrowbroad}Peak.fa",
        neg=GENOME,
        db=MEME_DB,
    output: 
        touch("results/{narrowbroad}_peaks_contrast_level/{contrast}/meme_clean_peaks/{contrast_name}.finished")
    log:
        "log/{narrowbroad}_peaks_contrast_level/{contrast}/meme_clean_peaks/{contrast_name}.log"
    benchmark:
        "log/{narrowbroad}_peaks_contrast_level/{contrast}/meme_clean_peaks/{contrast_name}.benchmark"
    params:
        odir="results/{narrowbroad}_peaks_contrast_level/{contrast}/meme_clean_peaks/",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    threads:
        6
    conda:
        "../envs/meme.yaml"
    shell:
        """
        meme-chip -oc {params.odir} -meme-p {threads} -db {input.db} {input.fasta} &> {log}
        """


## SUMMIT REGION ##

rule get_summit_neighbour:
# todo: filter blacklist, cpm filter
    input:
        summit="results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_summits.bed",
        genome=GENOME
    output:
        "results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_summits.{width}.fa"
    log:
        "log/narrow_peaks_contrast_level/{contrast}/{contrast_name}_summits.{width}.fa.log"
    threads:
        1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 6000
    conda:
        "../envs/biopython.yaml"
    shell:
        """
        python workflow/scripts/get_summit_neighbour.py {input.genome} {input.summit} {wildcards.width} {output} &> {log}
        """

rule meme_neibour:
    input: 
        fasta="results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_summits.{width}.fa",
        neg=GENOME,
        db=MEME_DB,
    output: 
        touch("results/narrow_peaks_contrast_level/{contrast}/memeSummit.{width}/{contrast_name}.finished")
    log:
        "log/narrow_peaks_contrast_level/{contrast}/memeSummit.{width}/{contrast_name}.log"
    benchmark:
        "log/narrow_peaks_contrast_level/{contrast}/memeSummit.{width}/{contrast_name}.benchmark"
    params:
        odir="results/narrow_peaks_contrast_level/{contrast}/memeSummit.{width}/",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    threads:
        6
    conda:
        "../envs/meme.yaml"
    shell:
        """
        meme-chip -oc {params.odir} -meme-p {threads} -db {input.db} {input.fasta} &> {log}
        """


## SUMMIT REGION BY CHR ##

rule split_fa_by_chr:
    input:
        fasta="results/narrow_peaks_contrast_level/{contrast}/{contrast_name}_summits.{width}.fa"
    output:
        "results/narrow_peaks_contrast_level/{contrast}/by_chr/{contrast_name}_summits.{width}.{chr}.fa"
    log:
        "log/narrow_peaks_contrast_level/{contrast}/by_chr/{contrast_name}_summits.{width}.{chr}.fa.log"
    threads:
        1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    shell:
        """
        grep -A 1 '>{wildcards.chr}:' {input} > {output} 2> {log}
        """


rule meme_neibour_chr_split:
    input: 
        fasta=rules.split_fa_by_chr.output,
        neg=GENOME,
        db=MEME_DB,
    output: 
        touch("results/narrow_peaks_contrast_level/{contrast}/memeSummit_chr.{width}_{chr}/{contrast_name}.finished")
    log:
        "log/narrow_peaks_contrast_level/{contrast}/memeSummit_chr.{width}_{chr}/{contrast_name}.log"
    benchmark:
        "log/narrow_peaks_contrast_level/{contrast}/memeSummit_chr.{width}_{chr}/{contrast_name}.benchmark"
    params:
        odir="results/narrow_peaks_contrast_level/{contrast}/memeSummit_chr.{width}_{chr}/",
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    threads:
        6
    conda:
        '../envs/meme.yaml'
    shell:
        """
        meme-chip -oc {params.odir} -meme-p {threads} -db {input.db} {input.fasta} &> {log}
        """
