import os


localrules:
    flair_combine_bed12,
    flair_transcriptomes

rule flair:
    input:
        expand("results/flair/transcriptome/{tissue}.gtf", tissue=util.tissues),


rule flair_bam_to_bed12:
    input:
        bam=input_long_read_bam,
        bai=input_long_read_bai,
    output:
        bed12="resources/bed12/{sample}.bed12",
    log:
        "logs/flair/bam_to_bed12/{sample}.log",
    threads: 1
    conda:
        "../envs/flair.yaml"
    shell:
        "(bam2Bed12 -i {input.bam} > {output.bed12}) > {log} 2>&1"


rule flair_combine_bed12:
    input:
        bed12=lambda wildcards: expand(
            "resources/bed12/{sample}.bed12",
            sample=util.samples_for_tissue(wildcards.tissue),
        ),
    output:
        "results/flair/{tissue}_combined.bed12",
    threads: 1
    conda:
        "../envs/flair.yaml"
    shell:
        "cat {input.bed12} > {output}"


rule flair_correct:
    input:
        bed12="results/flair/{tissue}_combined.bed12",
        gtf="resources/annotation.gtf",
        ref_fa="resources/reference.fa",
    output:
        correct="results/flair/correct/{tissue}_all_corrected.bed",
        inconsistent="results/flair/correct/{tissue}_all_inconsistent.bed",
        # "results/flair/correct/{tissue}_cannot_verify.bed" not always generated
    log:
        "logs/flair/correct/{tissue}.log",
    params:
        output_prefix=lambda wc, output: output["correct"].replace("_all_corrected.bed", ""),
    threads: 8
    resources:
        mem_mb=16 * 1024,
        runtime_min=60,
    conda:
        "../envs/flair.yaml"
    shell:
        "flair correct --query {input.bed12} --genome {input.ref_fa} --gtf {input.gtf} --threads {threads} --output {params.output_prefix} > {log} 2>&1"


rule flair_collapse:
    input:
        bed12="results/flair/correct/{tissue}_all_corrected.bed",
        gtf="resources/annotation.gtf",
        reads=lambda wildcards: util.longreads_for_tissue(wildcards.tissue),
        ref_fa="resources/reference.fa",
    output:
        # "results/flair/collapse/{tissue}.isoforms.bed",
        temp("results/flair/collapse/{tissue}.isoforms.gtf"),
        # "results/flair/collapse/{tissue}.isoforms.fa",
    log:
        "logs/flair/collapse/{tissue}_all_collapsed.log",
    params:
        output_prefix=lambda wc, output: output[0].replace(f".isoforms.gtf", ""),
    threads: 8
    resources:
        mem_mb=64 * 1024,
        runtime_min=24 * 60,
    conda:
        "../envs/flair.yaml"
    script:
        "../scripts/flair/collapse.py"


rule flair_transcriptomes:
    input:
        gtf="results/flair/collapse/{tissue}.isoforms.gtf",
    output:
        "results/flair/transcriptome/{tissue}.gtf",
    shell:
        "cp {input} {output}"
