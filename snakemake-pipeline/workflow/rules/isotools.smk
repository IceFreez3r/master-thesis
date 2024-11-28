rule isotools:
    input:
        expand("results/isotools_v0/transcriptome/{tissue}.gtf", tissue=util.tissues),


rule isotools_create:
    input:
        bams=lambda wildcards: util.long_read_bam_for_tissue(wildcards.tissue),
        bais=lambda wildcards: util.long_read_bai_for_tissue(wildcards.tissue),
        annotation_gff="resources/annotation_sorted.gtf.gz",
        annotation_tbi="resources/annotation_sorted.gtf.gz.tbi",
        ref_fa="resources/reference.fa",
        ref_fai="resources/reference.fa.fai",
        sample_table=config["sample_table"],
    output:
        gtf="results/isotools{version}/transcriptome/{tissue}.gtf",
        pkl="results/isotools{version}/pkl/{tissue}.pkl",
    log:
        "logs/isotools{version}/{tissue}.log",
    params:
        samples=lambda wildcards: util.samples_for_tissue(wildcards.tissue),
        extra=lambda wildcards: config["isotools"]["extra"].get('isotools' + wildcards.version, {}),
    threads: 1
    resources:
        mem_mb=128 * 1024,
        runtime_min=6 * 60,
    conda:
        # Uses advanced filters, which aren't available in the public pip version -> custom environment
        "isotools"
    script:
        "../scripts/isotools/create.py"
