rule CAGE:
    input:
        expand("resources/CAGE/{tissue}.bed", tissue=tissues)


rule unzip_CAGE_bed:
    output:
        "resources/CAGE/{tissue}.bed"
    params:
        dir = config["CAGE_dir"]
    shell:
        "gunzip -c {params.dir}/*{wildcards.tissue}*.bed.gz > {output}"
