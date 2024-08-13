rule isotools:
    input:
        expand("results/isotools/transcriptome/{tissue}.gtf", tissue=util.tissues),


rule isotools_create:
    input:
        bams=lambda wildcards: [input_long_read_bam({"sample": sample}) for sample in util.samples_for_tissue(wildcards.tissue)],
        bais=lambda wildcards: [input_long_read_bai({"sample": sample}) for sample in util.samples_for_tissue(wildcards.tissue)],
        annotation_gff=config["annot_gff"],
        annotation_tbi=config["annot_gff"] + ".tbi",
        reference_fa="resources/reference.fa",
        reference_fai="resources/reference.fa.fai",
        sample_table=config["sample_table"],
    output:
        gtf="results/isotools{conda}/transcriptome/{tissue}.gtf",
    log:
        "logs/isotools{conda}/{tissue}.log",
    params:
        samples=util.samples,
        query=config["isotools"]["query"],
    threads: 1
    resources:
        mem_mb=128 * 1024,
        runtime_min=6 * 60,
    conda:
        # Uses advanced filters, which aren't available in the public pip version -> custom environment
        "isotools{conda}"
    script:
        "../scripts/isotools/create.py"
