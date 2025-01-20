def isotools_versions():
    return [tool for tool in WORKING_TOOLS if "isotools" in tool]

def isotools_training_sequences_versions():
    # All versions that have "training sequences" defined
    return [tool for tool in isotools_versions() if config["isotools"].get(tool).get('training_sequences')]

rule isotools:
    input:
        expand("results/{tool}/transcriptome/{group}.gtf", tool=isotools_versions(), group=util.tissues),
        expand("results/{tool}/pkl/{group}.pkl", tool=isotools_versions(), group=util.tissues),
        expand("results/{tool}/training_sequences/{group}_{filter}.fa", tool=isotools_training_sequences_versions(),
               group=util.tissues, filter=['positive', 'negative']),


rule isotools_create:
    input:
        bams=lambda wildcards: util.long_read_bam_for_tissue(wildcards.group),
        bais=lambda wildcards: util.long_read_bai_for_tissue(wildcards.group),
        annotation_gff="resources/annotation_sorted.gtf.gz",
        annotation_tbi="resources/annotation_sorted.gtf.gz.tbi",
        ref_fa="resources/reference.fa",
        ref_fai="resources/reference.fa.fai",
        sample_table=config["sample_table"],
    output:
        gtf="results/isotools{version}/transcriptome/{group}.gtf",
        pkl="results/isotools{version}/pkl/{group}.pkl",
    log:
        "logs/isotools{version}/{group}.log",
    params:
        samples=lambda wildcards: util.samples_for_tissue(wildcards.group),
        extra=lambda wildcards: config["isotools"].get('isotools' + wildcards.version),
    threads: 1
    resources:
        mem_mb=128 * 1024,
        runtime_min=6 * 60,
    conda:
        # Uses advanced filters, which aren't available in the public pip version -> custom environment
        "isotools"
    script:
        "../scripts/isotools/create.py"


rule isotools_training_sequences:
    input:
        pkl="results/isotools{version}/pkl/{group}.pkl",
        ref_fa="resources/reference.fa",
        ref_fai="resources/reference.fa.fai",
        sqanti_annotation="results/sqanti/isotools{version}/qc/{group}/{group}_classification.txt",
    output:
        positive="results/isotools{version}/training_sequences/{group}_positive.fa",
        negative="results/isotools{version}/training_sequences/{group}_negative.fa",
    log:
        "logs/isotools{version}/{group}_training_sequences.log",
    params:
        extra=lambda wildcards: config["isotools"].get('isotools' + wildcards.version).get('training_sequences'),
        output_prefix=lambda wildcards, output: output.positive.replace("_positive.fa", ""),
    threads: 1
    resources:
        mem_mb=128 * 1024,
        runtime_min=6 * 60,
    conda:
        "isotools"
    script:
        "../scripts/isotools/training_sequences.py"
