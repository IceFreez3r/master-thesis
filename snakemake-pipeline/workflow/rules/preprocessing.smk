# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/minimap2/index.html
rule minimap_index:
    input:
        target="resources/reference.fa",
    output:
        "resources/reference_genome.mmi",
    log:
        "logs/common/minimap_index.log",
    threads: 32
    resources:
        mem_mb=32 * 1024,
        runtime_min=60,
    params:
        extra="-x splice:hq",
    wrapper:
        "v3.12.1/bio/minimap2/index"


rule minimap_align_longreads:
    input:
        query=lambda wildcards: util.longreads_for_sample(wildcards.sample),
        target="resources/reference_genome.mmi",
    output:
        bam="resources/mapped_reads/{sample}.bam",
        bam_sorted="resources/mapped_reads/{sample}_sorted.bam",
    log:
        "logs/common/minimap2_align/longreads/{sample}.log",
    params:
        sorting="coordinate",
        extra="-ax splice:hq -uf --MD",
    conda:
        "../envs/minimap2.yaml"
    threads: 8
    resources:
        mem_mb=32 * 1024,
        runtime_min=60,
    shell:
        """
        (
            minimap2 -t {threads} {params.extra} {input.target} {input.query} -o {output.bam}
            samtools sort -O BAM -@ {threads} -o {output.bam_sorted} {output.bam}
        ) > {log} 2>&1
        """


# rule index_gff_annotation:
#     input:
#         gff=config["annot_gff"],
#     output:
#         tbi=f"{config['annot_gff']}.tbi",
#     shell:
#         "tabix -p gff ${refdir}/${gff}_sorted.gff3.gz"
