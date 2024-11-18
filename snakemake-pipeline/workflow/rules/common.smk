import pandas as pd
import glob
import os

localrules:
    preprocess_reference_fa,
    unzip_rnaseq_reads,
    preprocess_polyA_peaks,
    samtools_index

WORKING_TOOLS = config["tools"]

if (config["test_run"]):
    for overwrite in config["test_config"]:
        config[overwrite] = config["test_config"][overwrite]


class Utility:
    def __init__(self, config):
        self.config = config
        self.sample_df = pd.read_csv(config['sample_table'], sep="\t")
        self.use_short_reads = config["rnaseq_fastq_fofn"] is not None
        self.cage_table = pd.read_csv(config["CAGE_table"], sep="\t")

        # Remove all samples from tissues without CAGE support
        self.sample_df = self.sample_df[self.sample_df["group"].isin(self.tissues)]
        if self.use_short_reads:
            self.rnaseq_fastq = pd.read_csv(config["rnaseq_fastq_fofn"], sep="\t")
            self.rnaseq_fastq = self.rnaseq_fastq[self.rnaseq_fastq["sample ID"].isin(self.samples)]
        else:
            self.rnaseq_fastq = pd.DataFrame(columns=["sample ID", "experiment ID", "file"])

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
        cage_tissues = self.cage_table["group"].unique().tolist()
        tissues = [tissue for tissue in tissues if tissue in cage_tissues]
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

    def experiment_for_sample(self, sample):
        experiment = self.sample_df[self.sample_df["sample ID"] == sample]["file"].values[0]
        experiment = experiment.split("/")[-1].split(".")[0]
        return experiment

    def tissue_for_sample(self, sample):
        tissue = self.sample_df[self.sample_df["sample ID"] == sample]["group"].values[0]
        return tissue

    def long_read_bam_for_sample(self, sample):
        sample_file = self.sample_df[self.sample_df["sample ID"] == sample]["file"].values[0]
        if sample_file.endswith(".bam"):
            return sample_file
        if config["test_run"]:
            return f"resources/test/mapped_reads/{wildcards["sample"]}_sorted_test.bam"
        return f"resources/mapped_reads/{wildcards["sample"]}_sorted.bam"

    def long_read_bai_for_sample(self, sample):
        sample_file = self.sample_df[self.sample_df["sample ID"] == sample]["file"].values[0]
        if sample_file.endswith(".bam"):
            return f"{sample_file}.bai"
        if config["test_run"]:
            return f"resources/test/mapped_reads/{wildcards["sample"]}_sorted_test.bam.bai"
        return f"resources/mapped_reads/{wildcards["sample"]}_sorted.bam.bai"


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

    def CAGE_file_for_tissue(self, tissue):
        return self.cage_table[self.cage_table["group"] == tissue]["file"].values[0]

    def long_read_bam_for_tissue(self, tissue):
        sample_files = self.sample_df[self.sample_df["group"] == tissue]["file"]
        if sample_files.iloc[0].endswith(".bam"):
            return sample_files.tolist()
        samples = self.samples_for_tissue(tissue)
        if config["test_run"]:
            return [f"resources/test/mapped_reads/{sample}_sorted_test.bam" for sample in samples]
        return [f"resources/mapped_reads/{sample}_sorted.bam" for sample in samples]


    def long_read_bai_for_tissue(self, tissue):
        sample_files = self.sample_df[self.sample_df["group"] == tissue]["file"]
        if sample_files.iloc[0].endswith(".bam"):
            return [f"{sample_file}.bai" for sample_file in sample_files]
        samples = self.samples_for_tissue(tissue)
        if config["test_run"]:
            return [f"resources/test/mapped_reads/{sample}_sorted_test.bam.bai" for sample in samples]
        return [f"resources/mapped_reads/{sample}_sorted.bam.bai" for sample in samples]



util = Utility(config)


class ENCODE_data:
    def __init__(self):
        self.fragment_sizes = self.import_fragment_sizes()
        self.has_error = False

    def import_fragment_sizes(self):
        if not os.path.exists("resources/fragment_sizes.tsv"):
            return pd.DataFrame(columns=["sample", "fragment_size", "fragment_sd"])
        return pd.read_csv("resources/fragment_sizes.tsv", sep="\t")

    def export_fragment_sizes(self):
        os.makedirs("resources", exist_ok=True)
        self.fragment_sizes.to_csv("resources/fragment_sizes.tsv", sep="\t", index=False)

    @property
    def samples(self):
        return self.fragment_sizes["sample"].tolist()

    def fragment_mean(self, sample):
        if not sample in self.samples:
            self.fragment_size_from_encode(sample)
        return self.fragment_sizes.loc[self.fragment_sizes["sample"] == sample, "fragment_size"].values[0]

    def fragment_sd(self, sample):
        if not sample in self.samples:
            self.fragment_size_from_encode(sample)
        return self.fragment_sizes.loc[self.fragment_sizes["sample"] == sample, "fragment_sd"].values[0]

    def fragment_size_from_encode(self, sample):
        if not self.has_error:
            url = f"https://www.encodeproject.org/experiments/{sample}/?format=json"
            try:
                response = requests.get(url)
            except requests.exceptions.SSLError as e:
                print(f"[WARN] Error getting fragment size from ENCODE for {sample}: {e}")
                print("[WARN] Defaulting to fragment size 200 and sd 30 for all samples")
                self.has_error = True
            try:
                response.raise_for_status()
                data = response.json()
                fragment_size = int(data["replicates"][0]["library"]["average_fragment_size"])
                fragment_sd = int(data["replicates"][0]["library"]["fragment_length_CV"])
                self.fragment_sizes = pd.concat([self.fragment_sizes, pd.DataFrame({"sample": [sample], "fragment_size": [fragment_size], "fragment_sd": [fragment_sd]})])
                self.export_fragment_sizes()
                return
            except requests.exceptions.HTTPError as e:
                print(f"[WARN] Error getting fragment size from ENCODE for {sample}: {e}")
                print("[WARN] Defaulting to fragment size 200 and sd 30 for all samples")
                self.has_error = True
        self.fragment_sizes = pd.concat([self.fragment_sizes, pd.DataFrame({"sample": [sample], "fragment_size": [200], "fragment_sd": [30]})])

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
    conda:
        "../envs/tabix.yaml"
    shell:
        "tabix -p gff {input} > {log} 2>&1"


def get_overlap_tools(wildcards):
    tools = WORKING_TOOLS
    # Only take first isotools tool
    first = [tool for tool in tools if "isotools" in tool][0]
    tools = [tool for tool in tools if not "isotools" in tool or tool == first]
    return tools

def tissue_gtfs(wildcards):
    tools = get_overlap_tools(wildcards)
    return {f"{tool}{ext}": f"results/{tool}/transcriptome/{wildcards.tissue}_sorted.gtf.gz{ext}" for tool in tools for ext in ["", ".tbi"]}

rule overlap:
    input:
        expand("results/plots/{tissue}/upset_all.png", tissue=util.tissues),
        expand("results/plots/{tissue}/upset_filtered.png", tissue=util.tissues),


rule tool_overlap:
    input:
        unpack(tissue_gtfs),
    output:
        upset = "results/plots/upset/all/{tissue}.png",
        upset_filtered = "results/plots/upset/filtered/{tissue}.png",
        counts = "resources/overlap/{tissue}.tsv",
    log:
        "logs/common/tool_overlap/{tissue}.log"
    params:
        tools = get_overlap_tools,
        tss_error = 20,
        pas_error = 20,
        junction_error = 5,
    threads: 32
    resources:
        runtime_min = 4 * 60,
        mem_mb = 32 * 1024,
    conda:
        "../envs/upset.yaml"
    script:
        "../scripts/common/tool_overlap.py"


rule tool_overlap_mean:
    input:
        expand("resources/overlap/{tissue}.tsv", tissue=util.tissues),
    output:
        upset = "results/plots/upset/all/mean.png",
        upset_filtered = "results/plots/upset/filtered/mean.png",
    log:
        "logs/common/tool_overlap_mean.log"
    conda:
        "../envs/upset.yaml"
    script:
        "../scripts/common/tool_overlap_mean.py"
