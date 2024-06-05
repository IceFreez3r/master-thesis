rule unzip_annotation:
    output:
        "resources/annotation.gtf",
    params:
        gtf_gz=config["annot_gtf"],
    shell:
        "gunzip -c {params.gtf_gz} > {output}"


rule unzip_transcriptome:
    output:
        gtfs=expand("resources/transcriptome/{sample}.gtf", sample=util.samples),
    params:
        samples=config["transcriptome_table"],
    run:
        df = pd.read_csv(params["samples"], sep="\t")
        for i, row in df.iterrows():
            sample = row["sample ID"]
            gtf_gz = row["gtf"]
            shell("gunzip -c {gtf_gz} > resources/transcriptome/{sample}.gtf")


rule unzip_reference_annotation:
    input:
        gtf_gz=config["annot_gtf"],
    output:
        "resources/reference_annotation.gtf",
    shell:
        "gunzip -c {input.gtf_gz} > {output}"


# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/minimap2/index.html
rule minimap_index:
    input:
        target=config["reference_fa"],
    output:
        "resources/reference_genome.mmi",
    log:
        "logs/common/minimap_index.log",
    threads: 32
    params:
        extra= "-x splice:hq"
    wrapper:
        "v3.11.0/bio/minimap2/index"


# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/minimap2/aligner.html
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
        extra="-ax splice:hq -uf --MD"
    conda:
        "../envs/minimap2.yaml"
    threads: 8
    resources:
        mem_mib=32 * 1024,
        runtime_min=60,
    shell:
        """
        (
            minimap2 -t {threads} {params.extra} {input.target} {input.query} -o {output.bam}
            samtools sort -O BAM -@ {threads} -o {output.bam_sorted} {output.bam}
        ) > {log} 2>&1
        """

rule index_bam:
    input:
        bam="resources/mapped_reads/{sample}_sorted.bam"
    output:
        bai="resources/mapped_reads/{sample}_sorted.bam.bai"
    shell:
        "samtools index {input.bam}"

# rule index_gff_annotation:
#     input:
#         gff=config["annot_gff"],
#     output:
#         tbi=f"{config['annot_gff']}.tbi",
#     shell:
#         "tabix -p gff ${refdir}/${gff}_sorted.gff3.gz"
