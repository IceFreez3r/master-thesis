import os


rule sqanti:
    input:
        expand("results/sqanti/{tool}/qc/{tissue}/{tissue}_SQANTI3_report.html", tool=WORKING_TOOLS, tissue=util.tissues)


rule sqanti_qc:
    input:
        gtf = "results/{tool}/transcriptome/{tissue}.gtf",
        ref_gtf = "resources/annotation.gtf",
        cage = "resources/CAGE/{tissue}.bed",
        ref_fa = "resources/reference.fa",
        polyA_motifs = config["sqanti"]["polyA_motif_list"],
        polyA_peaks = "resources/PolyASitePeaks.bed",
        kallisto_expression = lambda wildcards: expand("results/kallisto/{{tool}}/{sample}/abundance.tsv", sample=util.rnaseq_samples_for_tissue(wildcards.tissue)),
        # fastq=lambda wildcards: expand("resources/rnaseq/{sample}.fastq", sample=util.rnaseq_samples_for_tissue(wildcards.tissue)),
        # fastq_gz=lambda wildcards: util.rnaseq_reads_for_tissue(wildcards.tissue),
        # rnaseq_fastq = "results/sqanti/rnaseq_fofn/{tissue}.fofn",
        # rnaseq_fastq_gz = "results/sqanti/rnaseq_fofn/{tissue}_gz.fofn",
        sj_tabs = lambda wildcards: expand("results/star/{sample}/{sample}.SJ.out.tab", sample=util.rnaseq_samples_for_tissue(wildcards.tissue)),
        SR_bams = lambda wildcards: expand("results/star/{{tissue}}/{sample}.Aligned.sortedByCoord.out.bam", sample=util.rnaseq_samples_for_tissue(wildcards.tissue)),
    output:
        "results/sqanti/{tool}/qc/{tissue}/{tissue}_SQANTI3_report.html",
        "results/sqanti/{tool}/qc/{tissue}/{tissue}_SQANTI3_report.pdf",
        "results/sqanti/{tool}/qc/{tissue}/{tissue}_classification.txt",
        "results/sqanti/{tool}/qc/{tissue}/{tissue}_corrected.fasta",
        "results/sqanti/{tool}/qc/{tissue}/{tissue}_junctions.txt",
        "results/sqanti/{tool}/qc/{tissue}/{tissue}.params.txt",
        gtf="results/sqanti/{tool}/qc/{tissue}/{tissue}_corrected.gtf",
        # rnaseq_fofn = "results/sqanti/rnaseq_fofn/{tool}/{tissue}.fofn"
    log:
        out = "logs/sqanti/{tool}/qc/{tissue}.log",
        error = "logs/sqanti/{tool}/qc/{tissue}.error.log"
    params:
        sqanti_qc = os.path.join(config["sqanti"]["path"], "sqanti3_qc.py"),
        output_dir = lambda wc, output: os.path.dirname(output.gtf),
        extra = "--report both --skipORF",
        extra_user = config["sqanti"]["extra"],
        expression = lambda wildcards, input: ",".join(input.kallisto_expression),
        sj_tabs = lambda wildcards, input: ",".join(input.sj_tabs),
        SR_bam_dir = lambda wildcards, input: os.path.dirname(input.SR_bams[0]),
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
                --coverage {params.sj_tabs} --SR_bam {params.SR_bam_dir} --expression {params.expression}\
                --polyA_motif_list {input.polyA_motifs} --polyA_peak {input.polyA_peaks}\
                -d {params.output_dir}\
                {params.extra} {params.extra_user} --cpus {threads}\
                {input.gtf} {input.ref_gtf} {input.ref_fa}
            df -h $MXQ_JOB_TMPDIR
        ) > {log.out} 2> {log.error}
        """

def get_sqanti_tools(wildcards):
    if wildcards["version"] == "isotools":
        return [tool for tool in WORKING_TOOLS if "isotools" in tool]
    else:
        return [tool for tool in WORKING_TOOLS if not "isotools" in tool or "isotools" + wildcards["version"] == tool]

rule sqanti_comparison_plots:
    input:
        classifications = expand("results/sqanti/{tool}/qc/{tissue}/{tissue}_classification.txt", tool=get_sqanti_tools, tissue=util.tissues),
    output:
        "results/plots/sqanti/{version}/transcript_counts.png",
        "results/plots/sqanti/{version}/transcript_counts_subcategory.png",
        "results/plots/sqanti/{version}/transcript_counts_subcategory_ISM.png",
        "results/plots/sqanti/{version}/CAGE_support.png",
        "results/plots/sqanti/{version}/TSS_ratio.png",
        "results/plots/sqanti/{version}/PolyA_site.png",
        "results/plots/sqanti/{version}/PolyA_motif.png",
        "results/plots/sqanti/{version}/CAGE_support_FSM.png",
        "results/plots/sqanti/{version}/CAGE_support_ISM.png",
        "results/plots/sqanti/{version}/CAGE_support_NIC.png",
        "results/plots/sqanti/{version}/CAGE_support_NNC.png",
        "results/plots/sqanti/{version}/CAGE_support_non_FSM.png",
        "results/plots/sqanti/{version}/polyA_site_FSM.png",
        "results/plots/sqanti/{version}/polyA_site_ISM.png",
        "results/plots/sqanti/{version}/polyA_site_NIC.png",
        "results/plots/sqanti/{version}/polyA_site_NNC.png",
        "results/plots/sqanti/{version}/polyA_motif_FSM.png",
        "results/plots/sqanti/{version}/polyA_motif_ISM.png",
        "results/plots/sqanti/{version}/polyA_motif_NIC.png",
        "results/plots/sqanti/{version}/polyA_motif_NNC.png",
        "results/plots/sqanti/{version}/CAGE_support_monoexons.png",
        "results/plots/sqanti/{version}/CAGE_support_ISM_no_monoexons.png",
        "results/plots/sqanti/{version}/CAGE_support_no_monoexons.png",
        "results/plots/sqanti/{version}/CAGE_support_no_monoexons_no_3fragment.png",
        "results/plots/sqanti/{version}/CAGE_support_FSM_no_monoexons_no_3fragment.png",
        "results/plots/sqanti/{version}/CAGE_support_ISM_no_monoexons_no_3fragment.png",
        "results/plots/sqanti/{version}/CAGE_support_NIC_no_monoexons_no_3fragment.png",
        "results/plots/sqanti/{version}/CAGE_support_NNC_no_monoexons_no_3fragment.png",
        stats = "results/plots/sqanti/{version}/stats.tsv",
    log:
        "logs/sqanti/comparison/{version}.log",
    conda:
        "../envs/seaborn.yaml"
    params:
        tissues = util.tissues,
        tools = lambda wildcards: get_sqanti_tools(wildcards),
        toolnames = lambda wildcards: [tool.split("_")[0] for tool in get_sqanti_tools(wildcards)] if wildcards["version"] != "isotools" else get_sqanti_tools(wildcards),
        classifications = lambda wildcards: {
            tool: {
                tissue: f"results/sqanti/{tool}/qc/{tissue}/{tissue}_classification.txt"
                for tissue in util.tissues
            }
            for tool in get_sqanti_tools(wildcards)
        },
        output_dir = "results/plots/sqanti/{version}",
        plot_titles = config["sqanti"]["plot_titles"],
    script:
        "../scripts/sqanti/sqanti_comparison_plots.py"
