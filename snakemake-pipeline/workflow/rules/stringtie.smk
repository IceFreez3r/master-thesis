rule stringtie:
    input:
        expand("results/stringtie/transcriptome/{tissue}.gtf", tissue=util.tissues),


rule stringtie_run:
    input:
        bam=lambda wildcards: util.long_read_bam_for_sample(sample),
        # bais=expand("resources/mapped_reads/{sample}_sorted.bam.bai", sample=util.samples),
        annot_gff="resources/annotation.gtf",
    output:
        gtf="results/stringtie/samples/{sample}.gtf",
    log:
        "logs/stringtie/samples/{sample}.log",
    params:
        stringtie_exe=config["stringtie"]["path"],
        extra=config["stringtie"]["extra"]["run"],
    threads: 8
    shell:
        "{params.stringtie_exe} {input.bam} -o {output.gtf} -p {threads} -L -G {input.annot_gff} {params.extra} > {log} 2>&1"


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
        extra=config["stringtie"]["extra"]["merge"],
    threads: 16
    shell:
        "{params.stringtie_exe} --merge -o {output.gtf} -p {threads} -G {input.annot_gff} {input.sample_gtfs} {params.extra} > {log} 2>&1"


# Some entries of the output gtf had no strand info, so they get filtered out
rule stringtie_filter:
    input:
        gtf="results/stringtie/merged/{tissue}.gtf",
    output:
        gtf="results/stringtie/transcriptome/{tissue}.gtf",
    log:
        "logs/stringtie/filter/{tissue}.log",
    conda:
        "../envs/pyranges.yaml"
    resources:
        disk_mb = lambda wc, input: max(4 * input.size_mb, 1000)
    script:
        "../scripts/stringtie/filter.py"
