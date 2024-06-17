import pandas as pd
import glob
import os

class Utility:
    def __init__(self, config):
        self.config = config
        self.sample_df = pd.read_csv(config['sample_table'], sep="\t")
        self.rnaseq_df = pd.read_csv(config["rnaseq_table"], sep="\t")

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

    def rnaseq_alignment_for_sample(self, sample):
        alignment = self.rnaseq_df[self.rnaseq_df["sample ID"] == sample]["file"].values[0]
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

    def longreads_for_tissue(self, tissue):
        return self.sample_df[self.sample_df["group"] == tissue]["file"].tolist()

    def rnaseq_alignments_for_tissue(self, tissue):
        samples = self.samples_for_tissue(tissue)
        return self.rnaseq_df[self.rnaseq_df["sample ID"].isin(samples)]["file"].tolist()

util = Utility(config)

rule samtools_index:
    input:
        "{sample}.bam"
    output:
        "{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule tissue_gtfs:
    input:
        expand(
            "results/{tool}/tissue_gtfs.fofn",
            tool=["flair", "isotools"],
        ),
    output:
        "results/tissues_gtf.fofn",
    shell:
        "cat {input} > {output}"
