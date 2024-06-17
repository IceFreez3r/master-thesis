import os


localrules:
    flair_tissue_gtfs,


rule flair:
    input:
        expand("results/flair/gtf/{tissue}.gtf", tissue=util.tissues),
        # expand("results/flair/collapse/{tissue}.isoforms.bed", tissue=util.tissues),
        # expand("results/flair/collapse/{tissue}.isoforms.fa", tissue=util.tissues),


rule bam_to_bed12:
    input:
        bam="resources/mapped_reads/{sample}_sorted.bam",
    output:
        bed12="resources/bed12/{sample}.bed12",
    log:
        "logs/flair/align/bam_to_bed12_{sample}.log",
    threads: 1
    resources:
        mem_mib=32 * 1024,
        runtime_min=60,
    conda:
        "../envs/flair.yaml"
    shell:
        """
        echo "Aligning {input.bam} to bed12..." > {log}
        (bam2Bed12 -i {input.bam} > {output.bed12}) > {log} 2>&1
        """


rule combine_bed12:
    input:
        bed12=lambda wildcards: expand(
            "resources/bed12/{sample}.bed12",
            sample=util.samples_for_tissue(wildcards.tissue),
        ),
    output:
        "results/flair/{tissue}_combined.bed12",
    threads: 1
    resources:
        mem_mib=32 * 1024,
        runtime_min=60,
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
    threads: 8
    resources:
        mem_mib=16 * 1024,
        runtime_min=60,
    conda:
        "../envs/flair.yaml"
    shell:
        "flair correct --query {input.bed12} --genome {params.reference_fa} --gtf {input.gtf} --threads {threads} --output {params.output_prefix} > {log} 2>&1"


rule flair_collapse:
    input:
        bed12="results/flair/correct/{tissue}_all_corrected.bed",
        gtf="resources/reference_annotation.gtf",
        reads=lambda wildcards: util.longreads_for_tissue(wildcards.tissue),
    output:
        "results/flair/collapse/{tissue}.isoforms.bed",
        "results/flair/collapse/{tissue}.isoforms.gtf",
        "results/flair/collapse/{tissue}.isoforms.fa",
    log:
        "logs/flair/collapse/{tissue}_all_collapsed.log",
    params:
        reference_fa=config["reference_fa"],
        output_prefix=lambda wildcards: "results/flair/collapse/" + wildcards.tissue,
    threads: 8
    resources:
        mem_mib=8 * 1024,
        runtime_min=60,
    conda:
        "../envs/flair.yaml"
    shell:
        "flair collapse --query {input.bed12} --genome {params.reference_fa} --reads {input.reads} --gtf {input.gtf} --threads {threads} --output {params.output_prefix} > {log} 2>&1"


rule flair_tissue_gtfs:
    input:
        gtf=expand("results/flair/collapse/{tissue}.isoforms.gtf", tissue=util.tissues),
    output:
        "results/flair/tissue_gtfs.fofn",
    run:
        # tsv with tissue and file path
        with open(output[0], "w") as f:
            for tissue in util.tissues:
                path = os.path.join(os.getcwd(), "results/flair/collapse", f"{tissue}.isoforms.gtf")
                f.write(f"flair\t{tissue}\t{path}\n")
