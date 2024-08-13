import pandas as pd
import glob
import os

localrules:
    preprocess_reference_fa,
    unzip_rnaseq_reads,
    preprocess_polyA_peaks,
    samtools_index

WORKING_TOOLS = ["flair", "isoquant", "isotools", "isotools_new_tss", "stringtie"]

if (config["test_run"]):
    for overwrite in config["test_config"]:
        config[overwrite] = config["test_config"][overwrite]


class Utility:
    def __init__(self, config):
        self.config = config
        self.sample_df = pd.read_csv(config['sample_table'], sep="\t")
        self.rnaseq_fastq = pd.read_csv(config["rnaseq_fofn"]["fastq"], sep="\t")
        self.rnaseq_bam = pd.read_csv(config["rnaseq_fofn"]["bam"], sep="\t")

        # Remove all samples from tissues without CAGE support
        self.sample_df = self.sample_df[self.sample_df["group"].isin(self.tissues)]
        self.rnaseq_fastq = self.rnaseq_fastq[self.rnaseq_fastq["sample ID"].isin(self.samples)]
        self.rnaseq_bam = self.rnaseq_bam[self.rnaseq_bam["sample ID"].isin(self.samples)]

    @property
    def samples(self):
        samples = self.sample_df["sample ID"].tolist()
        return samples

    @property
    def rnaseq_samples(self):
        return self.rnaseq_fastq["sample ID"].tolist()

    @property
    def tissues(self):
        tissues = self.sample_df["group"].unique().tolist()
        CAGE_dir = self.config["CAGE_dir"]
        tissues = [tissue for tissue in tissues if glob.glob(f"{CAGE_dir}/*{tissue}*.bed.gz")]
        return tissues

    @property
    def experiments(self):
        experiments = self.sample_df["file"].tolist()
        experiments = [experiment.split("/")[-1].split(".")[0] for experiment in experiments]
        return experiments

    def longreads_for_sample(self, sample):
        reads = self.sample_df[self.sample_df["sample ID"] == sample]["file"].values[0]
        return reads

    def rnaseq_sample_for_sample(self, sample):
        return self.rnaseq_fastq[self.rnaseq_fastq["sample ID"] == sample]["experiment ID"].values[0]

    def rnaseq_reads_for_sample(self, sample):
        return self.rnaseq_fastq[self.rnaseq_fastq["sample ID"] == sample]["file"].values[0]

    def rnaseq_alignment_for_sample(self, sample):
        alignment = self.rnaseq_bam[self.rnaseq_bam["sample ID"] == sample]["file"].values[0]
        return alignment

    def experiment_for_sample(self, sample):
        experiment = self.sample_df[self.sample_df["sample ID"] == sample]["file"].values[0]
        experiment = experiment.split("/")[-1].split(".")[0]
        return experiment

    def tissue_for_sample(self, sample):
        tissue = self.sample_df[self.sample_df["sample ID"] == sample]["group"].values[0]
        return tissue


    def samples_for_tissue(self, tissue):
        samples = self.sample_df[self.sample_df["group"] == tissue]["sample ID"].tolist()
        return samples

    def rnaseq_samples_for_tissue(self, tissue):
        samples = self.samples_for_tissue(tissue)
        return self.rnaseq_fastq[self.rnaseq_fastq["sample ID"].isin(samples)]["sample ID"].tolist()

    def longreads_for_tissue(self, tissue):
        return self.sample_df[self.sample_df["group"] == tissue]["file"].tolist()

    def rnaseq_reads_for_tissue(self, tissue):
        samples = self.samples_for_tissue(tissue)
        return self.rnaseq_fastq[self.rnaseq_fastq["sample ID"].isin(samples)]["file"].tolist()

    def rnaseq_alignments_for_tissue(self, tissue):
        samples = self.samples_for_tissue(tissue)
        return self.rnaseq_bam[self.rnaseq_bam["sample ID"].isin(samples)]["file"].tolist()

util = Utility(config)


class ENCODE_data:
    def __init__(self):
        self.fragment_sizes = {}
        self.has_error = False

    def fragment_mean(self, sample):
        if not sample in self.fragment_sizes:
            self.fragment_size_from_encode(sample)
        return self.fragment_sizes[sample][0]

    def fragment_sd(self, sample):
        if not sample in self.fragment_sizes:
            self.fragment_size_from_encode(sample)
        return self.fragment_sizes[sample][1]

    def fragment_size_from_encode(self, sample):
        if not self.has_error:
            url = f"https://www.encodeproject.org/experiments/{sample}/?format=json"
            response = requests.get(url)
            try:
                response.raise_for_status()
                data = response.json()
                fragment_size = data["replicates"][0]["library"]["average_fragment_size"]
                fragment_sd = data["replicates"][0]["library"]["fragment_length_CV"]
                self.fragment_sizes[sample] = (fragment_size, fragment_sd)
                return
            except requests.exceptions.HTTPError:
                print(f"[WARN] Error getting fragment size from ENCODE for {sample}")
                print("[WARN] Defaulting to fragment size 200 and sd 30 for all samples")
                self.has_error = True
        self.fragment_sizes[sample] = (200, 30)

encode = ENCODE_data()

def shadow_large_rules():
    return "copy-minimal" if config["shadow_large_rules"] else None

rule preprocess_reference_fa:
    """Removes scaffolds from the reference"""
    input:
        config["reference_fa"],
    output:
        "resources/reference.fa",
    run:
        with open(input[0]) as f, open(output[0], "w") as out:
            copy = True
            for line in f:
                if line.startswith(">"):
                    copy = line.startswith(">chr")
                if copy:
                    out.write(line)


rule unzip_annotation:
    input:
        gz = config["annot_gtf"]
    output:
        gtf = "resources/annotation.gtf"
    log:
        "logs/common/unzip_annotation.log"
    resources:
        disk_mb = 5 * 1024
    shell:
        "(gunzip -c {input.gz} > {output.gtf}) > {log} 2>&1"


rule index_reference:
    input:
        "resources/reference.fa"
    output:
        "resources/reference.fa.fai"
    log:
        "logs/common/index_reference.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input} > {log} 2>&1"


rule samtools_index:
    input:
        "{sample}.bam"
    output:
        "{sample}.bam.bai"
    log:
        "logs/common/samtools_index/{sample}.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input} > {log} 2>&1"


rule unzip_rnaseq_reads:
    input:
        lambda wildcards: util.rnaseq_reads_for_sample(wildcards.sample)
    output:
        "resources/rnaseq/{sample}.fastq"
    log:
        "logs/common/unzip_rnaseq_reads/{sample}.log"
    shell:
        "(gunzip -c {input} > {output}) > {log} 2>&1"


rule preprocess_polyA_peaks:
    """unzips file and changes chromosomes from '1' to 'chr1'"""
    input:
        gz = config["PolyASitePeaks"]
    output:
        bed = "resources/PolyASitePeaks.bed"
    log:
        "logs/common/preprocess_polyA_peaks.log"
    shell:
        "(gunzip -c {input.gz} | sed 's/^/chr/' - > {output.bed}) > {log} 2>&1"

def input_long_read_bam(wildcards):
    if config["test_run"]:
        return f"resources/test/mapped_reads/{wildcards["sample"]}_sorted_test.bam"
    return f"resources/mapped_reads/{wildcards["sample"]}_sorted.bam"

def input_long_read_bai(wildcards):
    if config["test_run"]:
        return f"resources/test/mapped_reads/{wildcards["sample"]}_sorted_test.bam.bai"
    return f"resources/mapped_reads/{wildcards["sample"]}_sorted.bam.bai"

rule test_long_read_bam:
    input:
        bam = "resources/mapped_reads/{sample}_sorted.bam",
        bai = "resources/mapped_reads/{sample}_sorted.bam.bai",
    output:
        "resources/test/mapped_reads/{sample}_sorted_test.bam"
    log:
        "logs/common/test_long_read_bam/{sample}.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        "(samtools view -b {input.bam} chr15 > {output}) > {log} 2>&1"


rule sort_gzip_gtf:
    input:
        "results/{tool}/transcriptome/{tissue}.gtf"
    output:
        "results/{tool}/transcriptome/{tissue}_sorted.gtf.gz"
    log:
        "logs/common/sort_gzip_gtf/{tool}_{tissue}.log"
    shell:
        """
        (grep -v "^#" {input} | sort -k1,1 -k4,4n | bgzip -c > {output}) > {log} 2>&1
        """


rule index_gtf:
    input:
        "{path}_sorted.gtf.gz"
    output:
        "{path}_sorted.gtf.gz.tbi"
    log:
        "logs/common/index_gtf/{path}.log"
    shell:
        "tabix -p gff {input} > {log} 2>&1"
