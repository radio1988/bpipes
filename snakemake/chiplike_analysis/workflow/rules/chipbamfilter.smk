DATA_TYPE = config['DATA_TYPE']

# rule markDup:
#     # same speed as bwa_map, slow
#     input:
#         bam="results/sorted_reads/{sample}.bam",
#         bai="results/sorted_reads/{sample}.bam.bai"
#     output:
#         #temp(
#         bam=temp("results/markDup/{sample}.bam"),
#         metrics="results/markDup/{sample}.markDup_metrics.txt",
#         bai=temp("results/markDup/{sample}.bam.bai"),
#         #)  # make temp and save storage
#     log:
#         "log/markDup/{sample}.markDup.log"
#     benchmark:
#         "log/markDup/{sample}.benchmark"
#     conda:
#         "../envs/chiplike.yaml"
#     threads:
#         2
#     resources:
#         mem_mb=lambda wildcards, attempt: attempt * 16000
#     shell:
#         """
#         module load picard/2.17.8
#         PICARD=/share/pkg/picard/2.17.8/picard.jar
        
#         java -Xmx30g -XX:ParallelGCThreads=2 -jar $PICARD MarkDuplicates \
#         I={input.bam} \
#         O={output.bam} \
#         M={output.metrics} \
#         REMOVE_DUPLICATES=true \
#         ASSUME_SORTED=true \
#         &> {log}

#         samtools index {output.bam} &>> {log}
#         """

rule mark_duplicates:
    input:
        bam="results/sorted_reads/{sample}.bam",
        bai="results/sorted_reads/{sample}.bam.bai"
    # optional to specify a list of BAMs; this has the same effect
    # of marking duplicates on separate read groups for a sample
    # and then merging
    output:
        bam=temp("results/markDup/{sample}.bam"),
        metrics="results/markDup/{sample}.markDup_metrics.txt",
    log:
        "log/markDup/{sample}.markDup.log"
    benchmark:
        "log/markDup/{sample}.benchmark"
    params:
        "REMOVE_DUPLICATES=true ASSUME_SORTED=true",
    threads:
        2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000
    wrapper:
        "v0.75.0/bio/picard/markduplicates"


# todo: new MODE, DATA_TYPE pattern
if DATA_TYPE == 'DamID':
    rule DamID_filter:
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
            python workflow/scripts/filter_bam.py {input.bam} {GENOME} GATC {output} &> {log}
            """
elif DATA_TYPE == 'ATAC':
    rule ATAC_filter:
        input:
            "results/markDup/{sample}.bam"
        output:
            "results/clean_reads/{sample}.bam"
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

            samtools index {output} &> {log}
            """
elif DATA_TYPE == 'ChIP':
    rule ChIP_filter:
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
            1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 8000
        conda:
            "../envs/chiplike.yaml"
        shell:
            """
            cp {input.bam} {output.bam} &> {log}
            samtools index {input.bam}
            """
else: 
    sys.exit("DATA_TYPE error, see config.yaml for details")

