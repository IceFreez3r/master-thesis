import os


rule flair:
    input:
        expand("results/flair/collapse/{tissue}.isoforms.bed", tissue=util.tissues),
        expand("results/flair/collapse/{tissue}.isoforms.gtf", tissue=util.tissues),
        expand("results/flair/collapse/{tissue}.isoforms.fa", tissue=util.tissues),


rule bam_to_bed12:
    input:
        bam=lambda wildcards: util.get_alignment_for_sample(wildcards.sample),
    output:
        bed12="resources/bed12/{sample}.bed12",
    log:
        "logs/flair/align/bam_to_bed12_{sample}.log",
    conda:
        "../envs/flair.yaml"
    shell:
        "(bam2Bed12 -i {input.bam} > {output.bed12}) > {log} 2>&1"


rule combine_bed12:
    input:
        bed12=lambda wildcards: expand(
            "resources/bed12/{sample}.bed12",
            sample=util.get_samples_for_tissue(wildcards.tissue),
        ),
    output:
        "results/flair/{tissue}_combined.bed12",
    conda:
        "../envs/flair.yaml"
    shell:
        "cat {input.bed12} > {output}"


rule flair_correct:
    input:
        bed12="results/flair/{tissue}_combined.bed12",
        gtf="resources/reference_annotation.gtf",
    output:
        "results/flair/correct/{tissue}_all_corrected.bed",
        "results/flair/correct/{tissue}_all_inconsistent.bed",
        # "results/flair/correct/{tissue}_cannot_verify.bed" not always generated
    log:
        "logs/flair/correct/{tissue}_all_corrected.log",
    params:
        reference_fa=config["reference_fa"],
        output_prefix=lambda wildcards: "results/flair/correct/" + wildcards.tissue,
    conda:
        "../envs/flair.yaml"
    threads: 8
    shell:
        "flair correct --query {input.bed12} --genome {params.reference_fa} --gtf {input.gtf} --threads {threads} --output {params.output_prefix} > {log} 2>&1"


rule flair_collapse:
    input:
        bed12="results/flair/correct/{tissue}_all_corrected.bed",
        gtf="resources/reference_annotation.gtf",
        reads=lambda wildcards: expand(
            os.path.join(config["reads_dir"], "{experiment}.fastq.gz"),
            experiment=util.get_longreads_for_tissue(wildcards.tissue),
        ),
    output:
        "results/flair/collapse/{tissue}.isoforms.bed",
        "results/flair/collapse/{tissue}.isoforms.gtf",
        "results/flair/collapse/{tissue}.isoforms.fa",
    log:
        "logs/flair/collapse/{tissue}_all_collapsed.log",
    params:
        reference_fa=config["reference_fa"],
        output_prefix=lambda wildcards: "results/flair/collapse/" + wildcards.tissue,
    conda:
        "../envs/flair.yaml"
    threads: 8
    shell:
        "flair collapse --query {input.bed12} --genome {params.reference_fa} --reads {input.reads} --gtf {input.gtf} --threads {threads} --output {params.output_prefix} > {log} 2>&1"
