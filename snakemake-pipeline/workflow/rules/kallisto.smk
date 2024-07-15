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
        "logs/{tool}/transcriptome/{tissue}.log"
    conda:
        "../envs/gffread.yaml"
    shell:
        "gffread -w {output.fa} -g {input.ref_fa} {input.gtf}"


rule kallisto_index:
    input:
        fasta="results/{tool}/transcriptome/{tissue}.fa"
    output:
        index="results/kallisto/{tool}/{tissue}/reference.kidx"
    log:
        "logs/kallisto/{tool}/{tissue}/index.log"
    threads: 8
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto index -i {output.index} -t {threads} {input.fasta} > {log} 2>&1"


# tuples of mean and sd fragment size for each sample
fragment_sizes = {}

def fragment_mean(sample):
    if not sample in fragment_sizes:
        fragment_size_from_encode(sample)
    return fragment_sizes[sample][0]

def fragment_sd(sample):
    if not sample in fragment_sizes:
        fragment_size_from_encode(sample)
    return fragment_sizes[sample][1]

def fragment_size_from_encode(sample):
    url = f"https://www.encodeproject.org/experiments/{sample}/?format=json"
    response = requests.get(url)
    response.raise_for_status()
    data = response.json()
    fragment_sizes[sample] = (data["replicates"][0]["library"]["average_fragment_size"],
                              data["replicates"][0]["library"]["fragment_length_CV"])


rule kallisto_quant:
    input:
        index=lambda wildcards: f"results/kallisto/{wildcards.tool}/{util.tissue_for_sample(wildcards.sample)}/reference.kidx",
        fastq=lambda wildcards: util.rnaseq_reads_for_sample(wildcards.sample),
    output:
        abundance="results/kallisto/{tool}/{sample}/abundance.tsv",
        h5="results/kallisto/{tool}/{sample}/abundance.h5",
        json="results/kallisto/{tool}/{sample}/run_info.json",
    log:
        "logs/kallisto/{tool}/{sample}.log"
    params:
        output_dir="results/kallisto/{tool}/{sample}",
        fragment_length=lambda wildcards: fragment_mean(util.rnaseq_sample_for_sample(wildcards.sample)),
        sd=lambda wildcards: fragment_sd(util.rnaseq_sample_for_sample(wildcards.sample)),
    threads: 8
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto quant -i {input.index} -o {params.output_dir} -t {threads} --single -l {params.fragment_length} -s {params.sd} {input.fastq} > {log} 2>&1"
