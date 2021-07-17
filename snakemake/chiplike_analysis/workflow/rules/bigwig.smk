BW_BIN_SIZE=config['BW_BIN_SIZE']
MODE=config['MODE']
minFragmentLength=config['minFragmentLength']
maxFragmentLength=config['maxFragmentLength']

def bamCoverage_pse_params(config):
    if MODE == 'PE':
        return "--minFragmentLength " + str(minFragmentLength) + " --maxFragmentLength " + str(maxFragmentLength) + " -e 150"
    elif MODE == 'SE':
        return ""
    else:
        print("MODE not PE nor SE:", MODE)

rule bamCoverage:
    # for ChIP
    input:
        "results/clean_reads/{sample}.bam"
    output:
        "results/clean_reads_bigWig/{sample}.cpm.bw"
    params:
        bin_size=BW_BIN_SIZE,
        pse=bamCoverage_pse_params
    threads:
        8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    log:
        "log/clean_reads_bigWig/{sample}.bamCoverage.log"
    benchmark:
        "log/clean_reads_bigWig/{sample}.bamCoverage.benchmark"
    conda:
        "../envs/deeptools.yaml"
    shell:
        # Aim: same as our downstream filters, extensions
        """
        bamCoverage --bam {input} \
        -o  {output} \
        --numberOfProcessors {threads} \
        --outFileFormat bigwig \
        --normalizeUsing CPM \
        --binSize {params.bin_size} \
        {params.pse} \
        &> {log}
        """


 
