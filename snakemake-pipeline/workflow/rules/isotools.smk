localrules:
    isotools_tissue_gtfs,


rule isotools:
    input:
        expand("results/isotools/gtf/{tissue}.gtf", tissue=util.tissues),


rule isotools_create:
    input:
        bams=expand("resources/mapped_reads/{sample}_sorted.bam", sample=util.samples),
        bais=expand(
            "resources/mapped_reads/{sample}_sorted.bam.bai", sample=util.samples
        ),
        annotation_gff=config["annot_gff"],
        annotation_tbi=config["annot_gff"] + ".tbi",
        reference_fa=config["reference_fa"],
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
        mem_mib=64 * 1024,
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
        tissue_gtfs=expand("results/isotools/gtf/{tissue}.gtf", tissue=util.tissues),
    log:
        "logs/isotools/tissues.log",
    params:
        tissues=util.tissues,
        query=config["isotools"]["query"],
        output_prefix="results/isotools/gtf/",
    threads: 1
    resources:
        mem_mib=16 * 1024,
        runtime_min=2 * 60,
    conda:
        # Uses advanced filters, which aren't available in the public pip version -> custom environment
        "isotools"
    script:
        "../scripts/isotools/tissues.py"


rule isotools_tissue_gtfs:
    input:
        gtf=expand("results/isotools/gtf/{tissue}.gtf", tissue=util.tissues),
    output:
        "results/isotools/tissue_gtfs.fofn",
    run:
        with open(output[0], "w") as f:
            for tissue in util.tissues:
                path = os.path.join(
                    os.getcwd(), "results/isotools/gtf", f"{tissue}.gtf"
                )
                f.write(f"isotools\t{tissue}\t{path}\n")
