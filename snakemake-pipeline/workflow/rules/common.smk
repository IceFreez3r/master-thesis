import pandas as pd
import glob
import os

localrules:
    preprocess_reference_fa,
    unzip_rnaseq_reads,
    preprocess_polyA_peaks,

WORKING_TOOLS = ["flair", "isotools", "isoquant", "stringtie"]

class Utility:
    def __init__(self, config):
        self.config = config
        self.sample_df = pd.read_csv(config['sample_table'], sep="\t")
        self.rnaseq_fastq = pd.read_csv(config["rnaseq_fofn"]["fastq"], sep="\t")
        self.rnaseq_bam = pd.read_csv(config["rnaseq_fofn"]["bam"], sep="\t")

    @property
    def samples(self):
        samples = self.sample_df["sample ID"].tolist()
        return samples

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

    def rnaseq_reads_for_sample(self, sample):
        reads = self.rnaseq_fastq[self.rnaseq_fastq["sample ID"] == sample]["file"].values[0]
        return reads

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
    shell:
        "(gunzip -c {input.gz} > {output.gtf}) > {log} 2>&1"


rule index_reference:
    input:
        "resources/reference.fa"
    output:
        "resources/reference.fa.fai"
    log:
        "logs/common/index_reference.log"
    shell:
        "samtools faidx {input} > {log} 2>&1"


rule samtools_index:
    input:
        "{sample}.bam"
    output:
        "{sample}.bam.bai"
    log:
        "logs/common/samtools_index/{sample}.log"
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
    '''unzips file and changes chromosomes from '1' to 'chr1''''
    input:
        gz = config["PolyASitePeaks"]
    output:
        bed = "resources/PolyASitePeaks.bed"
    log:
        "logs/common/preprocess_polyA_peaks.log"
    shell:
        "(gunzip -c {input.gz} | sed 's/^/chr/' - > {output.bed}) > {log} 2>&1"
