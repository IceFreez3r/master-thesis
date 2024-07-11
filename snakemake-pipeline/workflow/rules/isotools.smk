rule isotools:
    input:
        expand("results/isotools/transcriptome/{tissue}.gtf", tissue=util.tissues),


rule isotools_create:
    input:
        bams=expand("resources/mapped_reads/{sample}_sorted.bam", sample=util.samples),
        bais=expand(
            "resources/mapped_reads/{sample}_sorted.bam.bai", sample=util.samples
        ),
        annotation_gff=config["annot_gff"],
        annotation_tbi=config["annot_gff"] + ".tbi",
        reference_fa="resources/reference.fa",
        reference_fai="resources/reference.fa.fai",
        sample_table=config["sample_table"],
    output:
        pkl="results/isotools/isotools.pkl",
        gtf="results/isotools/isotools.gtf",
    log:
        "logs/isotools/create.log",
    params:
        samples=util.samples,
        tissues=util.tissues,
        coverage_threshold=config["isotools"]["coverage_threshold"],
    threads: 1
    resources:
        mem_mb=128 * 1024,
        runtime_min=6 * 60,
    conda:
        # Uses advanced filters, which aren't available in the public pip version -> custom environment
        "isotools"
    script:
        "../scripts/isotools/create.py"


rule isotools_tissues:
    input:
        pkl="results/isotools/isotools.pkl",
    output:
        tissue_gtfs=expand("results/isotools/transcriptome/{tissue}.gtf", tissue=util.tissues),
    log:
        "logs/isotools/tissues.log",
    params:
        tissues=util.tissues,
        query=config["isotools"]["query"],
        output_prefix="results/isotools/transcriptome/",
    threads: 1
    resources:
        mem_mb=128 * 1024,
        runtime_min=2 * 60,
    conda:
        # Uses advanced filters, which aren't available in the public pip version -> custom environment
        "isotools"
    script:
        "../scripts/isotools/tissues.py"
