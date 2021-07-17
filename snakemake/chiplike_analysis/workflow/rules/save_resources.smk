BLACKLIST=config['BLACKLIST']
GTF=config['GTF']

# to simplify threads/params/resources
localrules: save_blacklist, save_gtf

rule save_blacklist:
    input:
        BLACKLIST
    output:
        "resources/"+os.path.split(BLACKLIST)[1]+".gz"
    params:
        "resources/"+os.path.split(BLACKLIST)[1]
    shell:
        """
        cp {input} {params}
        gzip {params}
        """


rule save_gtf:
    input:
        GTF
    output:
        "resources/"+os.path.split(GTF)[1]+".gz"
    params:
        "resources/"+os.path.split(GTF)[1]
    shell:
        """
        cp {input} {params}
        gzip {params}
        """
