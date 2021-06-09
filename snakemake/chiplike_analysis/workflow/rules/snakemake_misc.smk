rule create_dag:
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000 
    threads:
        1
    output:
        "Workflow_DAG.svg"
    log:
        "log/create_dag/Workflow_DAG.svg.log"
    shell:
        "snakemake --dag all | dot -Tsvg > {output} 2> {log}"


rule reset:
    shell:
        """
        snakemake --unlock
        """