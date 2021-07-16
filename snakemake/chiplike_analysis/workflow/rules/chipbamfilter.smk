DATA_TYPE = config['DATA_TYPE']

rule markDup:
    # same speed as bwa_map, slow
    input:
        bam="results/sorted_reads/{sample}.bam",
        bai="results/sorted_reads/{sample}.bam.bai"
    output:
        bam=temp("results/markDup/{sample}.bam"),
        metrics="results/markDup/{sample}.markDup_metrics.txt",
    log:
        "log/markDup/{sample}.markDup.log"
    benchmark:
        "log/markDup/{sample}.benchmark"
    conda:
        "../envs/picard.yaml"
    threads:
        1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000

    shell:
        """
        # todo: portability, not only picard from conda (2.25.5) can be used
        p=`which picard` && echo $p &> {log};
        PICARD=`echo $p|sed 's/bin\/picard/share\/picard-2.25.5-0\/picard.jar/'` && echo $PICARD &>>{log}; 
        java -Xmx15g -jar $PICARD \
        MarkDuplicates \
        I={input.bam} \
        O={output.bam} \
        M={output.metrics} \
        REMOVE_DUPLICATES=true \
        ASSUME_SORTED=true \
        &>> {log}
        """


# todo: new MODE, DATA_TYPE pattern
if DATA_TYPE == 'DamID':
    rule bam_filter_damid:
        input:
            bam="results/sorted_reads/{sample}.bam",
            bai="results/sorted_reads/{sample}.bam.bai"
        output:
            "results/clean_reads/{sample}.bam"
        log:
            "log/clean_reads/{sample}.log"
        benchmark:
            "log/clean_reads/{sample}.benchmark"
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 16000
        threads:
            1
        conda:
            "../envs/chiplike.yaml"
        shell:
            """
            # need samtools/1.9
            python scripts/filter_bam.py {input.bam} {GENOME} GATC {output} &> {log}
            """
elif DATA_TYPE == 'ATAC':
    rule bam_filter_atac:
        input:
            "results/markDup/{sample}.bam"
        output:
            bam="results/clean_reads/{sample}.bam",
            bai="results/clean_reads/{sample}.bam.bai"
        log:
            "log/clean_reads/{sample}.log"
        benchmark:
            "log/clean_reads/{sample}.benchmark"
        threads:
            2
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 8000
        conda:
            "../envs/chiplike.yaml"
        shell:
            """
                echo 'ATACseq mode'
                echo {input}
                echo 'remove mito reads; keep paired reads with MQ>20 and 38-2000nt fragment size only'
                samtools view -h {input} 2>{log}| perl -lane 'print unless ($F[2] eq {chrM} and $_ != /\@/)' 2>>{log}| awk \'{config[filter]}\' 2>>{log}| $samtools sort -m 8G -o {output}  2>> {log}
                cp {input} {output}

                samtools index {output.bam} {output.bai} &>> {log}

            """
elif DATA_TYPE == 'ChIP':
    rule bam_filter_chip:
        input:
            bam="results/markDup/{sample}.bam",
        output:
            bam="results/clean_reads/{sample}.bam",
            bai="results/clean_reads/{sample}.bam.bai"
        log:
            "log/clean_reads/{sample}.log"
        benchmark:
            "log/clean_reads/{sample}.benchmark"
        threads:
            2
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 8000
        conda:
            "../envs/chiplike.yaml"
        shell:
            """
            cp {input.bam} {output.bam} &> {log}
            samtools index {output.bam} {output.bai} &>> {log}
            """
else: 
    sys.exit("DATA_TYPE error, see config.yaml for details")

