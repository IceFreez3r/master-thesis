import os


localrules:
    flair_combine_bed12,
    flair_transcriptomes,
    flair_fix

rule flair:
    input:
        expand("results/flair/transcriptome/{tissue}.gtf", tissue=util.tissues),


rule flair_bam_to_bed12:
    input:
        bam=lambda wildcards: util.long_read_bam_for_sample(sample),
        bai=lambda wildcards: util.long_read_bai_for_sample(sample),
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
        extra=config["flair"]["extra"]["correct"],
    threads: 8
    resources:
        mem_mb=16 * 1024,
        runtime_min=60,
    conda:
        "../envs/flair.yaml"
    shell:
        """
        flair correct --query {input.bed12} --genome {input.ref_fa} --gtf {input.gtf}\
              --threads {threads} --output {params.output_prefix} {params.extra} > {log} 2>&1
        """


rule flair_collapse:
    input:
        bed12="results/flair/correct/{tissue}_all_corrected.bed",
        gtf="resources/annotation.gtf",
        reads=lambda wildcards: util.longreads_for_tissue(wildcards.tissue),
        ref_fa="resources/reference.fa",
    output:
        gtf = temp("results/flair/collapse/{tissue}.isoforms.gtf"),
    log:
        "logs/flair/collapse/{tissue}_all_collapsed.log",
    params:
        output_prefix=lambda wc, output: output.gtf.replace(f".isoforms.gtf", ""),
        bed_split_size=config["flair"]["bed_split_size"],
        extra=config["flair"]["extra"]["collapse"],
    threads: 16
    resources:
        mem_mb=64 * 1024,
        runtime_min=24 * 60,
    conda:
        "../envs/flair.yaml"
    script:
        "../scripts/flair/collapse.py"


# Fix 0-length/negative length exons, required for indexing
rule flair_fix:
    input:
        gtf = "results/flair/collapse/{tissue}.isoforms.gtf",
    output:
        gtf = temp("results/flair/fixed/{tissue}.isoforms.gtf"),
    log:
        "logs/flair/fix/{tissue}.log",
    run:
        with open(input.gtf) as f, open(output.gtf, "w") as out, open(log[0], "w") as log:
            for line in f:
                if line.startswith("#"):
                    out.write(line)
                else:
                    fields = line.split("\t")
                    # FLAIR produces 0-length/negative length exons for lung
                    # See [GitHub issue](https://github.com/BrooksLabUCSC/flair/issues/356)
                    start = int(fields[3])
                    end = int(fields[4])
                    if start >= end:
                        log.write(f"Fixing {end - start} length exon: {line}")
                        fields[4] = str(start + 1)
                    out.write("\t".join(fields))


rule flair_transcriptomes:
    input:
        gtf="results/flair/fixed/{tissue}.isoforms.gtf",
    output:
        "results/flair/transcriptome/{tissue}.gtf",
    shell:
        "cp {input} {output}"
