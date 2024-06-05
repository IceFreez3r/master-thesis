import pandas as pd
import os
import sys


rule talon:
    input:
        expand("results/talon/talon_{tissue}.gtf", tissue=util.tissues),


rule talon_label_reads:
    input:
        bam=os.path.join(config["alignment_dir"], "{experiment}_aligned.bam"),
    output:
        sam="results/talon/labeled/{experiment}_labeled.sam",
        tsv="results/talon/labeled/{experiment}_read_labels.tsv",
    params:
        genome=config["reference_fa"],
        output_prefix="results/talon/labeled/{experiment}",
    log:
        "logs/talon/label_reads/{experiment}.log",
    threads: 8
    resources:
        mem_mib=32 * 1024,
        runtime_min=60,
    shadow:
        "shallow"  # call generates temp files
    conda:
        "../envs/talon.yaml"
    shell:
        "talon_label_reads --g {params.genome} --f {input.bam} --t {threads} --o {params.output_prefix} > {log} 2>&1"


rule talon_initialize_database:
    input:
        annot_gtf="results/annotation.gtf",
    output:
        "results/talon/talon.db",
    params:
        genome_name=config["genome_name"],
        annot_name=config["annot_name"],
        output_prefix="results/talon/talon",
    log:
        "logs/talon/init_db.log",
    threads: 8
    resources:
        mem_mib=32 * 1024,
        runtime_min=60,
    conda:
        "../envs/talon.yaml"
    shell:
        "talon_initialize_database --f {input.annot_gtf} --g {params.genome_name} --a {params.annot_name} --o {params.output_prefix} > {log} 2>&1"


# rule talon_config:
#     '''creates a csv file with dataset name, sample description, platform, sam file'''
#     input:
#         sample_sams = expand('results/talon/labeled/{experiment}_labeled.sam', experiment=util.experiments)
#     output:
#         csv = 'results/talon/config.csv'
#     params:
#         sequencing_platform = config['sequencing_platform']
#     run:
#         config = []
#         for sample in samples:
#             sample_index = samples[samples == sample].index[0]
#             config.append({
#                 'dataset_name': sample,
#                 'sample_description': sample_df.iloc[sample_index]['Description'],
#                 'platform': params['sequencing_platform'],
#                 'sam': input['sample_sams'][sample_index]
#             })
#         config = pd.DataFrame(config)
#         config.to_csv(output['csv'], index=False, header=False)


rule talon_run:
    input:
        db="results/talon/talon.db",
        sample_sams=expand(
            "results/talon/labeled/{experiment}_labeled.sam",
            experiment=util.experiments,
        ),
    output:
        config="results/talon/config.csv",
        log="results/talon/talon_QC.log",
        tsv="results/talon/read_annot.tsv",
    log:
        "logs/talon/talon.log",
    params:
        samples=util.samples,
        genome_name=config["genome_name"],
        output_prefix="results/talon/",
        use_tmpdir=config["talon"]["use_tmpdir"],
        sequencing_platform=config["sequencing_platform"],
    threads: 64
    resources:
        mem_mib=800 * 1024,
        runtime_min=lambda wildcards, attempt: 24 * 60 * attempt**2,
        disk_mib=800 * 1024,
    conda:
        "../envs/talon.yaml"
    script:
        "../scripts/talon.py"


rule talon_abundance:
    input:
        db="results/talon/talon.db",
        log="results/talon/talon_QC.log",  # ensures talon has finished
    output:
        "results/talon/abundance.tsv",
    log:
        "logs/talon/abundance.log",
    params:
        genome_name=config["genome_name"],
        annot_name=config["annot_name"],
        output_prefix="results/talon/",
    threads: 32
    resources:
        mem_mib=400 * 1024,
        runtime_min=lambda wildcards, attempt: 120 * attempt**2,
        disk_mib=400 * 1024,
    conda:
        "../envs/talon.yaml"
    shell:
        "talon_abundance --db {input.db} --annot {params.annot_name} --build {params.genome_name} --o {params.output_prefix} > {log} 2>&1"


rule talon_filter_transcripts:
    input:
        db="results/talon/talon.db",
        log="results/talon/talon_QC.log",  # ensures talon has finished
    output:
        tsv="results/talon/filtered_talon.tsv",
    log:
        "logs/talon/filter.log",
    params:
        annot_name=config["annot_name"],
        maxFracA=0.5,  # TODO: Move to config
        minCount=5,  # TODO: Move to config
        minDatasets=2,
    threads: 32
    resources:
        mem_mib=400 * 1024,
        runtime_min=lambda wildcards, attempt: 120 * attempt**2,
        disk_mib=400 * 1024,
    conda:
        "../envs/talon.yaml"
    shell:
        "talon_filter_transcripts --db {input.db} --annot {params.annot_name} --maxFracA {params.maxFracA} --minCount {params.minCount} --minDatasets {params.minDatasets} --o {output.tsv} > {log} 2>&1"


rule talon_create_GTF:
    input:
        db="results/talon/talon.db",
        log="results/talon/talon_QC.log",  # ensures talon has finished
    output:
        dataset="results/talon/talon_{tissue}.txt",
        gtf="results/talon/talon_{tissue}.gtf",
    log:
        "logs/talon/create_gtf_{tissue}.log",
    params:
        genome_name=config["genome_name"],
        annot_name=config["annot_name"],
        output_prefix=lambda wildcards: "results/talon/talon_{wildcards.tissue}",
        samples=lambda wildcards: "\n".join(util.samples_for_tissue(wildcards.tissue)),
    threads: 32
    resources:
        mem_mib=400 * 1024,
        runtime_min=lambda wildcards, attempt: 120 * attempt**2,
        disk_mib=400 * 1024,
    conda:
        "../envs/talon.yaml"
    shell:
        """
        echo {params.samples} > {output.dataset}
        talon_create_GTF --db {input.db} -annot {params.annot_name} --build {params.genome_name} --datasets {output.dataset} --o {params.output_prefix} > {log} 2>&1
        """


rule talon_create_adata:
    input:
        db="results/talon/talon.db",
        log="results/talon/talon_QC.log",  # ensures talon has finished
    output:
        h5ad="results/talon/talon.h5ad",
    log:
        "logs/talon/create_adata.log",
    params:
        genome_name=config["genome_name"],
        annot_name=config["annot_name"],
    threads: 32
    resources:
        mem_mib=400 * 1024,
        runtime_min=lambda wildcards, attempt: 120 * attempt**2,
        disk_mib=400 * 1024,
    conda:
        "../envs/talon.yaml"
    shell:
        "talon_create_adata --db {input.db} -annot {params.annot_name} --build {params.genome_name} --o {output.h5ad} > {log} 2>&1"
