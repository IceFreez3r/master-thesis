import requests


rule kallisto:
    input:
        abundance=expand("results/kallisto/{sample}/abundance.tsv", sample=util.rnaseq_samples_for_tissue),


rule kallisto_index:
    input:
        fasta=config["annot_fa"]
    output:
        index="results/kallisto/reference.kidx"
    log:
        "logs/kallisto/index.log"
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
        index="results/kallisto/reference.kidx",
        fastq=lambda wildcards: util.rnaseq_reads_for_sample(wildcards.sample),
    output:
        abundance="results/kallisto/{sample}/abundance.tsv",
        h5="results/kallisto/{sample}/abundance.h5",
        json="results/kallisto/{sample}/run_info.json",
    log:
        "logs/kallisto/{sample}.log"
    params:
        output_dir=lambda wildcards: f"results/kallisto/{wildcards.sample}",
        fragment_length=config["kallisto"]["fragment-length"],
        sd=config["kallisto"]["sd"]
    threads: 8
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto quant -i {input.index} -o {params.output_dir} -t {threads} --single -l {params.fragment_length} -s {params.sd} {input.fastq} > {log} 2>&1"
