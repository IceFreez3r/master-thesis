rule isotools:
    input:
        expand("results/isotools/{tissue}.gtf", tissue=util.tissues)

rule isotools_create:
    input:
        bams=expand("resources/mapped_reads/{sample}_sorted.bam", sample=util.samples),
        bais=expand("resources/mapped_reads/{sample}_sorted.bam.bai", sample=util.samples),
    output:
        pkl="resources/isotools/isotools.pkl",
        gtf="resources/isotools/isotools.gtf",
    log:
        "logs/isotools/create.log",
    params:
        reference_fa=config["reference_fa"],
        annotation_gff=config["annot_gff"],
        samples=util.samples,
        tissues=util.tissues,
        sample_table=config["sample_table"],
        coverage_threshold=config['isotools']["coverage_threshold"],
    conda:
        # Uses advanced filters, which aren't available in the public pip version -> custom environment
        "isotools",
    script:
        "../scripts/isotools/create.py"

rule isotools_tissues:
    input:
        pkl="resources/isotools/isotools.pkl",
    output:
        tissue_gtfs=expand("results/isotools/{tissue}.gtf", tissue=util.tissues),
    log:
        "logs/isotools/tissues.log",
    params:
        tissues=util.tissues,
        query=config['isotools']["query"],
    conda:
        # Uses advanced filters, which aren't available in the public pip version -> custom environment
        "isotools",
    script:
        "../scripts/isotools/tissues.py"
