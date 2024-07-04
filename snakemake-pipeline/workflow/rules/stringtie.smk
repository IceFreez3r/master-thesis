rule stringtie:
    input:
        expand("results/stringtie/tissues/{tissue}.gtf", tissue=util.tissues),


rule stringtie_run:
    input:
        bam="resources/mapped_reads/{sample}_sorted.bam",
        # bais=expand("resources/mapped_reads/{sample}_sorted.bam.bai", sample=util.samples),
        annot_gff="resources/reference_annotation.gtf",
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
        annot_gff="resources/reference_annotation.gtf",
    output:
        gtf="results/stringtie/transcriptome/{tissue}.gtf",
    log:
        "logs/stringtie/transcriptome/{tissue}.log",
    params:
        stringtie_exe=config["stringtie"]["path"],
    threads: 16
    shell:
        "{params.stringtie_exe} --merge -o {output.gtf} -p {threads} -G {input.annot_gff} {input.sample_gtfs} > {log} 2>&1"
