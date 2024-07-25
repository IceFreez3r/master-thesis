import requests
import os


rule kallisto:
    input:
        abundance=expand("results/kallisto/{tool}/{sample}/abundance.tsv", tool=WORKING_TOOLS, sample=util.rnaseq_samples),


rule tissue_gtf_to_fa:
    input:
        gtf="results/{tool}/transcriptome/{tissue}.gtf",
        ref_fa="resources/reference.fa",
        ref_fai="resources/reference.fa.fai",
    output:
        fa="results/{tool}/transcriptome/{tissue}.fa"
    log:
        "logs/kallisto/{tool}/gtf_to_fa/{tissue}.log"
    conda:
        "../envs/gffread.yaml"
    shell:
        "gffread -w {output.fa} -g {input.ref_fa} {input.gtf} > {log} 2>&1"


rule kallisto_index:
    input:
        fasta="results/{tool}/transcriptome/{tissue}.fa"
    output:
        index="results/kallisto/{tool}/index/{tissue}.kidx"
    log:
        "logs/kallisto/{tool}/index/{tissue}.log"
    threads: 16
    resources:
        mem_mb=16 * 1024,
        disk_mb=lambda wc, input: 10 * input.size_mb
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto index -i {output.index} -t {threads} {input.fasta} > {log} 2>&1"


rule kallisto_quant:
    input:
        index=lambda wildcards: f"results/kallisto/{wildcards.tool}/index/{util.tissue_for_sample(wildcards.sample)}.kidx",
        fastq=lambda wildcards: util.rnaseq_reads_for_sample(wildcards.sample),
    output:
        abundance="results/kallisto/{tool}/{sample}/abundance.tsv",
        h5="results/kallisto/{tool}/{sample}/abundance.h5",
        json="results/kallisto/{tool}/{sample}/run_info.json",
    log:
        "logs/kallisto/{tool}/{sample}.log"
    params:
        output_dir=lambda wc, output: output.abundance.replace("abundance.tsv", ""),
        fragment_length=lambda wildcards: encode.fragment_mean(util.rnaseq_sample_for_sample(wildcards.sample)),
        sd=lambda wildcards: encode.fragment_sd(util.rnaseq_sample_for_sample(wildcards.sample)),
    threads: 8
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto quant -i {input.index} -o {params.output_dir} -t {threads} --single -l {params.fragment_length} -s {params.sd} {input.fastq} > {log} 2>&1"
