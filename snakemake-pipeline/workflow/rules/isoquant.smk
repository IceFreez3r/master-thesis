import os


rule isoquant:
    input:
        expand("results/isoquant/{tissue}/OUT/OUT.transcript_models.gtf", tissue=util.tissues),


def check_for_annotation_db(wildcards):
    # isoquant preprocesses the annotation file
    # speed up the process by reusing the created .db file
    filename = os.path.basename(config["annot_gtf"]).split(".gz")[0]
    db_name = f"results/isoquant/{filename}.db"
    if os.path.exists(db_name):
        return db_name
    else:
        return config["annot_gtf"]


rule isoquant_run:
    input:
        bams=lambda wildcards: expand("resources/mapped_reads/{sample}_sorted.bam", sample=util.samples_for_tissue(wildcards.tissue)),
        bais=lambda wildcards: expand("resources/mapped_reads/{sample}_sorted.bam.bai", sample=util.samples_for_tissue(wildcards.tissue)),
        ref_fa="resources/reference.fa",
        annot_gtf=check_for_annotation_db,
    output:
        "results/isoquant/{tissue}/OUT/OUT.transcript_models.gtf"
    log:
        "logs/isoquant/{tissue}.log",
    params:
        output_folder=lambda wildcards: f"results/isoquant/{wildcards.tissue}",
    threads: 32
    resources:
        mem_mib=512 * 1024,
        runtime_min=8 * 60,
    conda:
        "../envs/isoquant.yaml"
    shell:
        "isoquant.py --reference {input.ref_fa} --genedb {input.annot_gtf} --bam {input.bams} --data_type pacbio_ccs -o {params.output_folder} --threads {threads} --complete_genedb --sqanti_output > {log} 2>&1"


rule isoquant_tissue_gtfs:
    input:
        expand("results/isoquant/{tissue}/OUT/OUT.transcript_models.gtf", tissue=util.tissues),
    output:
        "results/isoquant/tissue_gtfs.fofn",
    run:
        with open(output[0], "w") as f:
            for tissue in util.tissues:
                path = os.path.join(
                    os.getcwd(), f"results/isoquant/{tissue}", f"OUT/OUT.transcript_models.gtf"
                )
                f.write(f"isoquant\t{tissue}\t{path}\n")
