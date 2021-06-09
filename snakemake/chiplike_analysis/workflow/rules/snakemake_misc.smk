rule create_dag:
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000 
    threads:
        1
    output:
        "Workflow_DAG.svg"
    params:
        target_rule="targets"
    log:
        "log/create_dag/Workflow_DAG.svg.log"
    shell:
        "snakemake --dag {params.target_rule} | dot -Tsvg > {output} 2> {log}"


rule reset:
    shell:
        """
        snakemake --unlock
        """