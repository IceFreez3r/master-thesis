rule CAGE:
    input:
        expand("resources/CAGE/{tissue}.bed", tissue=util.tissues)


rule unzip_CAGE_bed:
    input:
        lambda wildcards: util.CAGE_file_for_tissue(wildcards.tissue)
    output:
        "resources/CAGE/{tissue}.bed"
    log:
        "logs/common/unzip_CAGE_bed/{tissue}.log"
    resources:
        disk_mb=lambda wildcards, input: input.size * 5
    shell:
        "(gunzip -c {input} > {output}) > {log} 2>&1"
