localrules:
    sqanti_rnaseq_fofn

rule sqanti:
    input:
        expand("results/sqanti/{tool}/qc/{tissue}/{tissue}_SQANTI3_report.html", tool=WORKING_TOOLS, tissue=util.tissues)

rule sqanti_rnaseq_fofn:
    input:
        fastq=lambda wildcards: expand("resources/rnaseq/{sample}.fastq", sample=util.rnaseq_samples_for_tissue(wildcards.tissue)),
    output:
        "results/sqanti/rnaseq_fofn/{tissue}.fofn"
    run:
        with open(output[0], 'w') as f:
            f.write("\n".join(input.fastq))

rule sqanti_qc:
    input:
        gtf="results/{tool}/transcriptome/{tissue}.gtf",
        ref_gtf = "resources/annotation.gtf",
        cage = "resources/CAGE/{tissue}.bed",
        ref_fa = "resources/reference.fa",
        polyA_motifs = config["sqanti"]["polyA_motif_list"],
        polyA_peaks = "resources/PolyASitePeaks.bed",
        kallisto_expression = lambda wildcards: expand("results/kallisto/{sample}/abundance.tsv", sample=util.rnaseq_samples_for_tissue(wildcards.tissue)),
        fastq=lambda wildcards: expand("resources/rnaseq/{sample}.fastq", sample=util.rnaseq_samples_for_tissue(wildcards.tissue)),
        rnaseq_fastq = "results/sqanti/rnaseq_fofn/{tissue}.fofn",
    output:
        "results/sqanti/{tool}/qc/{tissue}/{tissue}_SQANTI3_report.html",
        "results/sqanti/{tool}/qc/{tissue}/{tissue}_SQANTI3_report.pdf",
    log:
        out = "logs/sqanti/{tool}/qc/{tissue}.log",
        error = "logs/sqanti/{tool}/qc/{tissue}.error.log"
    params:
        sqanti_qc=os.path.join(config["sqanti"]["path"], "sqanti3_qc.py"),
        output_dir = "results/sqanti/{tool}/qc/{tissue}/",
        extra = "--report both",
        extra_user = config["sqanti"]["extra"],
        expression = lambda wildcards, input: ",".join(input.kallisto_expression),
    threads: 32
    resources:
        mem_mb=128 * 1024,
        runtime_min=24 * 60,
    conda:
        "../envs/sqanti.yaml"
    shell:
        """
        python {params.sqanti_qc} --CAGE_peak {input.cage}\
            --short_reads {input.rnaseq_fastq} --expression {params.expression}\
            --polyA_motif_list {input.polyA_motifs} --polyA_peak {input.polyA_peaks}\
            -d {params.output_dir}\
            {params.extra} {params.extra_user} --cpus {threads}\
            {input.gtf} {input.ref_gtf} {input.ref_fa} > {log.out} 2> {log.error}
        """

rule sqanti_filter:
    input:
        sqanti_class = "results/sqanti/qc/{tissue}/{tissue}.classification.txt", # TODO
    output:
        classification = "results/sqanti/filter/{filter_type}/{tissue}/{tissue}.classification.filtered.txt",
    params:
        output_dir = "results/sqanti/filter/{tissue}/",
        output_prefix = "results/sqanti/filter/{tissue}/"
        """
        TODO: Do I need these?
        --isoAnnotGFF3 ISOANNOTGFF3: Path to the isoAnnotLite GFF3 file that is to be filtered.
        --isoforms ISOFORMS: Path to the fasta/fastq isoform file that is to be filtered.
        --gtf GTF: Path to the GTF file that is to be filtered.
        --sam SAM: Path to the SAM alignment of the input fasta/fastq.
        --faa FAA: Path to the ORF prediction faa file to be filtered by SQANTI3.
        """
        # TODO: custom rule file?
    conda:
        "../envs/sqanti.yaml"
    shell:
        "python sqanti3_filter.py {wildcards.filter_type} {input.sqanti_class} -d {params.output_dir} -o {params.output_prefix}"
