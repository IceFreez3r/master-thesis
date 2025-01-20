import os


rule sqanti:
    input:
        expand("results/sqanti/{tool}/qc/{group}/{group}_SQANTI3_report.html", tool=WORKING_TOOLS, group=util.tissues)


def sqanti_short_read_input(wildcards):
    if not util.use_short_reads:
        return {}
    return {
        "sj_tabs": expand("results/star/{sample}/{sample}.SJ.out.tab", sample=util.rnaseq_samples_for_tissue(wildcards.group)),
        "sr_bams": expand("results/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam", sample=util.rnaseq_samples_for_tissue(wildcards.group)),
        "kallisto": expand("results/kallisto/{tool}/{sample}/abundance.tsv", tool=WORKING_TOOLS, sample=util.rnaseq_samples_for_tissue(wildcards.group)),
    }

rule sqanti_qc:
    input:
        unpack(sqanti_short_read_input),
        gtf = "results/{tool}/transcriptome/{group}.gtf",
        ref_gtf = "resources/annotation.gtf",
        cage = "resources/CAGE/{group}.bed",
        ref_fa = "resources/reference.fa",
        polyA_motifs = config["sqanti"]["polyA_motif_list"],
        polyA_peaks = "resources/PolyASitePeaks.bed",
    output:
        "results/sqanti/{tool}/qc/{group}/{group}_SQANTI3_report.html",
        "results/sqanti/{tool}/qc/{group}/{group}_SQANTI3_report.pdf",
        "results/sqanti/{tool}/qc/{group}/{group}_classification.txt",
        "results/sqanti/{tool}/qc/{group}/{group}_corrected.fasta",
        "results/sqanti/{tool}/qc/{group}/{group}_junctions.txt",
        "results/sqanti/{tool}/qc/{group}/{group}.params.txt",
        gtf="results/sqanti/{tool}/qc/{group}/{group}_corrected.gtf",
    log:
        out = "logs/sqanti/{tool}/qc/{group}.log",
        error = "logs/sqanti/{tool}/qc/{group}.error.log",
    params:
        sqanti_qc = os.path.join(config["sqanti"]["path"], "sqanti3_qc.py"),
        output_dir = lambda wc, output: os.path.dirname(output.gtf),
        extra = "--report both --skipORF",
        extra_user = config["sqanti"].get("extra", ""),
        short_read_params = lambda wildcards, input: f"--coverage {','.join(input.sj_tabs)} --SR_bam {os.path.dirname(input.sr_bams[0])} --expression {','.join(input.kallisto)}" if util.use_short_reads else "",
        # expression = lambda wildcards, input: ','.join(input.kallisto),
        # sj_tabs = lambda wildcards, input: ','.join(input.sj_tabs),
        # SR_bam_dir = lambda wildcards, input: os.path.dirname(input.sr_bams[0]),
    threads: 32
    resources:
        mem_mb=128 * 1024,
        runtime_min=24 * 60,
        # FLAIR transcriptomes tend to be a lot larger
        disk_mb=lambda wc, input: max(10*input.size_mb * (5 if wc.tool == "flair" else 1), 1000)
    conda:
        "../envs/sqanti.yaml"
    shell:
        """
        (
            python {params.sqanti_qc} --CAGE_peak {input.cage}\
                {params.short_read_params}\
                --polyA_motif_list {input.polyA_motifs} --polyA_peak {input.polyA_peaks}\
                -d {params.output_dir}\
                {params.extra} {params.extra_user} --cpus {threads}\
                {input.gtf} {input.ref_gtf} {input.ref_fa}
        ) > {log.out} 2> {log.error}
        """

class SQANTI_plots:
    def __init__(self):
        self.plot_groups = config["sqanti"]["plot_groups"]
        # Replace "all" with all tools
        for version, groups in self.plot_groups.items():
            if groups["tools"] == "all":
                self.plot_groups[version]["tools"] = WORKING_TOOLS
            if not "tool_names" in groups:
                self.plot_groups[version]["tool_names"] = groups["tools"]
            else:
                self.plot_groups[version]["tool_names"] = groups["tool_names"]
            assert len(self.plot_groups[version]["tools"]) == len(self.plot_groups[version]["tool_names"]), f"Number of tools and toolnames must match for {version}"
            if groups["groups"] == "all":
                self.plot_groups[version]["groups"] = util.tissues
            if not "group_names" in groups:
                self.plot_groups[version]["group_names"] = groups["groups"]
            else:
                self.plot_groups[version]["group_names"] = groups["group_names"]
            assert len(self.plot_groups[version]["groups"]) == len(self.plot_groups[version]["group_names"]), f"Number of groups and group_names must match for {version}"

    def __getitem__(self, version):
        return self.plot_groups[version]


sqanti_plots = SQANTI_plots()


rule sqanti_comparison_plots:
    input:
        classifications = lambda wildcards: expand("results/sqanti/{tool}/qc/{group}/{group}_classification.txt",
                                                   tool=sqanti_plots[wildcards.plot_group]["tools"],
                                                   group=sqanti_plots[wildcards.plot_group]["groups"]),
    output:
        "results/plots/sqanti/{plot_group}/transcript_counts.png",
        "results/plots/sqanti/{plot_group}/transcript_counts_subcategory.png",
        "results/plots/sqanti/{plot_group}/transcript_counts_subcategory_ISM.png",
        expand("results/plots/sqanti/{{plot_group}}/{TSS_metric}{transcript_category}{filter}.png",
               TSS_metric=["CAGE_support", "TSS_ratio"],
               transcript_category=["", "_FSM", "_ISM", "_NIC", "_NNC"],
               filter=["", "_no_3prime", "_no_monoexons", "_no_monoexons_no_3prime"]),
        expand("results/plots/sqanti/{{plot_group}}/{PAS_metric}{transcript_category}{filter}.png",
               PAS_metric=["polyA_site", "polyA_motif"],
               transcript_category=["", "_FSM", "_ISM", "_NIC", "_NNC"],
               filter=["", "_no_5prime", "_no_monoexons", "_no_monoexons_no_5prime"]),
        stats = "results/plots/sqanti/{plot_group}/stats.tsv",
    log:
        "logs/sqanti/comparison/{plot_group}.log",
    conda:
        "../envs/seaborn.yaml"
    params:
        tools = lambda wildcards: sqanti_plots[wildcards.plot_group]["tools"],
        tool_names = lambda wildcards: sqanti_plots[wildcards.plot_group]["tool_names"],
        groups = lambda wildcards: sqanti_plots[wildcards.plot_group]["groups"],
        group_names = lambda wildcards: sqanti_plots[wildcards.plot_group]["group_names"],
        classifications = lambda wildcards: {
            tool: {
                group: f"results/sqanti/{tool}/qc/{group}/{group}_classification.txt"
                for group in sqanti_plots[wildcards.plot_group]["groups"]
            }
            for tool in sqanti_plots[wildcards.plot_group]["tools"]
        },
        output_dir = lambda wildcards, output: os.path.dirname(output[0]),
        plot_titles = config["sqanti"]["plot_titles"],
        dpi = config["sqanti"]["dpi"],
        tss_cmap = config["sqanti"]["tss_cmap"],
        pas_cmap = config["sqanti"]["pas_cmap"],
    script:
        "../scripts/sqanti/sqanti_comparison_plots.py"
