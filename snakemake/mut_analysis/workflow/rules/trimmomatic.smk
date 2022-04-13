MODE=config['MODE']

if MODE == "PE":
    rule trimmomatic:
        input:
            r1="fastq/{sample}.R1.fastq.gz",
            r2="fastq/{sample}.R2.fastq.gz"
        output:
            r1="results/trimmed_reads/{sample}.R1.fastq.gz",
            r2="results/trimmed_reads/{sample}.R2.fastq.gz",
            u1="results/trimmed_reads/{sample}.R1.unpaired.fastq.gz",
            u2="results/trimmed_reads/{sample}.R2.unpaired.fastq.gz"
            # todo: SE
        params:
            trim='ILLUMINACLIP:resources/adapters/TruSeq3-PE.fa:2:30:10  \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 2000 
            # human need 18G
        threads:
            4
        log:
            "log/trimmomatic/{sample}.trimomatic.log"
        benchmark:
            "log/trimmomatic/{sample}.trimomatic.benchmark"
        conda:
            "../envs/trimmomatic.yaml"
        shell:
            """
            trimmomatic \
            PE -phred33 \
            {input.r1} {input.r2} \
            {output.r1} {output.u1} \
            {output.r2} {output.u2} \
            {params.trim}
            """
else:
    raise Error ("SE trimmomatic not implemented yet")
