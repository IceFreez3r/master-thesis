import os


localrules:
    isoquant_transcriptomes


rule isoquant:
    input:
        expand("results/isoquant/transcriptome/{tissue}.gtf", tissue=util.tissues),


def check_for_annotation_db(wildcards):
    # isoquant preprocesses the annotation file
    # speed up the process by reusing the created .db file
    filename = os.path.basename(config["annot_gtf"]).split(".gz")[0]
    db_name = os.path.join("results/isoquant", wildcards.tissue, "{filename}.db")
    if os.path.exists(db_name):
        return db_name
    else:
        return config["annot_gtf"]


rule isoquant_run:
    input:
        bams=lambda wildcards: util.long_read_bam_for_tissue(wildcards.tissue),
        bais=lambda wildcards: util.long_read_bai_for_tissue(wildcards.tissue),
        ref_fa="resources/reference.fa",
        annot_gtf=check_for_annotation_db,
    output:
        gtf=temp("results/isoquant/{tissue}/OUT/OUT.transcript_models.gtf"),
        db=f"results/isoquant/{{tissue}}/{os.path.basename(config["annot_gtf"]).split(".gz")[0]}.db",
    log:
        "logs/isoquant/{tissue}.log",
        "results/isoquant/{tissue}/isoquant.log",
    params:
        output_folder=lambda wc, output: output["gtf"].replace("/OUT/OUT.transcript_models.gtf", ""),
    threads: 32
    resources:
        mem_mb=512 * 1024,
        runtime_min=12 * 60,
        disk_mb=lambda wc, input: 4*input.size_mb,
    conda:
        "../envs/isoquant.yaml"
    shell:
        """
        isoquant.py --reference {input.ref_fa} --genedb {input.annot_gtf} --bam {input.bams} \
                    --data_type pacbio_ccs -o {params.output_folder} --threads {threads} \
                    --complete_genedb --sqanti_output > {log[0]} 2>&1
        """


rule isoquant_transcriptomes:
    input:
        gtf="results/isoquant/{tissue}/OUT/OUT.transcript_models.gtf",
    output:
        "results/isoquant/transcriptome/{tissue}.gtf",
    shell:
        "cp {input} {output}"
