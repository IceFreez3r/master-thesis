localrules:
    sqanti_rnaseq_fofn

rule sqanti:
    input:
        expand("results/sqanti/qc/this_wont_exist_{tool}_{tissue}.txt", tool=WORKING_TOOLS, tissue=util.tissues)

rule sqanti_rnaseq_fofn:
    input:
        fastq=lambda wildcards: expand("resources/rnaseq/{sample}.fastq", sample=util.rnaseq_sammples_for_tissue(wildcards.tissue)),
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
        fastq=lambda wildcards: expand("resources/rnaseq/{sample}.fastq", sample=util.rnaseq_sammples_for_tissue(wildcards.tissue)),
        rnaseq_fastq = "results/sqanti/rnaseq_fofn/{tissue}.fofn",
    output:
        "results/sqanti/qc/{tool}/{tissue}/{tissue}_SQANTI3_report.html"
    log:
        out = "logs/sqanti/qc/{tool}/{tissue}.log",
        error = "logs/sqanti/qc/{tool}/{tissue}.error.log"
    params:
        sqanti_qc=os.path.join(config["sqanti"]["path"], "sqanti3_qc.py"),
        # TODO: --fl full-length abundance data (from talon?)
        # TODO: One of the two is probably wrong/redundant
        output_dir = "results/sqanti/qc/{tool}/{tissue}/",
        # output_prefix = "results/sqanti/qc/{tool}/{tissue}/",
        extra = "--aligner_choice minimap2 --report both"
        extra_user = config["sqanti"]["extra"]
        # TODO: --force_id_ignore might be needed
    threads: 16
    resources:
        mem_mb=128 * 1024,
        runtime_min=8 * 60,
    conda:
        "../envs/sqanti.yaml"
    shell:
        """
        python {params.sqanti_qc} --CAGE_peak {input.cage}\
            --short_reads {input.rnaseq_fastq}\
            --polyA_motif_list {input.polyA_motifs}\
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
