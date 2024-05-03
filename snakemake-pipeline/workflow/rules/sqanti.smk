rule sqanti3_qc:
    input:
        gtf = "results/talon/gtfs/{sample}.gtf",
        ref_gtf = "results/annotation.gtf",
        # cage = "resources/CAGE/{tissue}.bed", # Include in analysis?
    # output:
    # TODO
    params:
        ref_fa = config["genome_fa"],
        # TODO: --CAGE_peak CAGE data
        # TODO: --short_reads short read data
        # TODO: --polyA_motif_list, provided by SQANTI repo?
        # TODO: --fl full-length abundance data (from talon?)
        output_dir = "results/sqanti3/qc/{sample}/",
        output_prefix = "results/sqanti3/qc/{sample}/",
        report_format = "both"
    threads: 8
    shell:
        # TODO: wrap in env?
        "python sqanti3_qc.py {input.gtf} {input.ref_gtf} {params.ref_fa} -d {params.output_dir} -o {params.output_prefix} --report {params.report_format} --cpus {threads}"

rule sqanti3_filter:
    input:
        sqanti_class = "results/sqanti3/qc/{sample}/{sample}.classification.txt", # TODO
    output:
        classification = "results/sqanti3/filter/{filter_type}/{sample}/{sample}.classification.filtered.txt",
    params:
        output_dir = "results/sqanti3/filter/{sample}/",
        output_prefix = "results/sqanti3/filter/{sample}/"
        """
        TODO: Do I need these?
        --isoAnnotGFF3 ISOANNOTGFF3: Path to the isoAnnotLite GFF3 file that is to be filtered.
        --isoforms ISOFORMS: Path to the fasta/fastq isoform file that is to be filtered.
        --gtf GTF: Path to the GTF file that is to be filtered.
        --sam SAM: Path to the SAM alignment of the input fasta/fastq.
        --faa FAA: Path to the ORF prediction faa file to be filtered by SQANTI3.
        """
        # TODO: custom rule file?
    shell:
        "python sqanti3_filter.py {wildcards.filter_type} {input.sqanti_class} -d {params.output_dir} -o {params.output_prefix}"
