rule stringtie:
    input:
        expand("results/stringtie/transcriptome/{tissue}.gtf", tissue=util.tissues),


rule stringtie_run:
    input:
        bam=input_long_read_bam,
        # bais=expand("resources/mapped_reads/{sample}_sorted.bam.bai", sample=util.samples),
        annot_gff="resources/annotation.gtf",
    output:
        gtf="results/stringtie/samples/{sample}.gtf",
    log:
        "logs/stringtie/samples/{sample}.log",
    params:
        stringtie_exe=config["stringtie"]["path"],
    threads: 8
    shell:
        "{params.stringtie_exe} {input.bam} -o {output.gtf} -p {threads} -L -G {input.annot_gff} > {log} 2>&1"


rule stringtie_merge:
    input:
        sample_gtfs=lambda wildcards: expand(
            "results/stringtie/samples/{sample}.gtf",
            sample=util.samples_for_tissue(wildcards.tissue),
        ),
        annot_gff="resources/annotation.gtf",
    output:
        gtf="results/stringtie/merged/{tissue}.gtf",
    log:
        "logs/stringtie/merged/{tissue}.log",
    params:
        stringtie_exe=config["stringtie"]["path"],
    threads: 16
    shell:
        "{params.stringtie_exe} --merge -o {output.gtf} -p {threads} -G {input.annot_gff} {input.sample_gtfs} > {log} 2>&1"


rule stringtie_filter:
    input:
        gtf="results/stringtie/merged/{tissue}.gtf",
    output:
        gtf="results/stringtie/transcriptome/{tissue}.gtf",
    log:
        "logs/stringtie/filter/{tissue}.log",
    conda:
        "../envs/pyranges.yaml"
    script:
        "../scripts/stringtie/filter.py"
