import pandas as pd
import os
import sys

sample_df = pd.read_csv(config['sample_table'], sep="\t")
samples = sample_df['file'].apply(lambda x: x.split("/")[-1].split('.')[0])

rule talon_label_reads:
    input:
        bam = os.path.join(config['alignment_dir'], '{sample}_aligned.bam'),
    output:
        sam = 'results/talon/labeled/{sample}_labeled.sam',
        tsv = 'results/talon/labeled/{sample}_read_labels.tsv'
    params:
        genome = config['genome_fa'],
        output_prefix = 'results/talon/labeled/{sample}'
    log:
        'logs/talon/label_reads/{sample}.log'
    threads: 8
    shadow: 'shallow' # call generates temp files
    conda: '../envs/talon.yaml'
    shell:
        'talon_label_reads --g {params.genome} --f {input.bam} --t {threads} --o {params.output_prefix} > {log} 2>&1'

rule talon_initialize_database:
    input:
        annot_gtf = 'results/annotation.gtf'
    output:
        'results/talon/talon.db'
    params:
        genome_name = config['genome_name'],
        annot_name = config['annot_name'],
        output_prefix = 'results/talon/talon'
    log:
        'logs/talon/init_db.log'
    conda: "../envs/talon.yaml"
    shell:
        'talon_initialize_database --f {input.annot_gtf} --g {params.genome_name} --a {params.annot_name} --o {params.output_prefix} > {log} 2>&1'

# rule talon_config:
#     '''creates a csv file with dataset name, sample description, platform, sam file'''
#     input:
#         sample_sams = expand('results/talon/labeled/{sample}_labeled.sam', sample=samples)
#     output:
#         csv = 'results/talon/config.csv'
#     params:
#         annot_gtf = config['annot_gtf'],
#         genome_name = config['genome_name'],
#         annot_name = config['annot_name'],
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

rule talon:
    input:
        db = 'results/talon/talon.db',
        sample_sams = expand('results/talon/labeled/{sample}_labeled.sam', sample=samples),
    output:
        config = 'results/talon/config.csv',
        log = 'results/talon/talon_QC.log'
    log:
        'logs/talon/talon.log'
    params:
        genome_name = config['genome_name'],
        output_prefix = 'results/talon/',
        use_tmpdir = config['talon_use_tmpdir']
    threads: 32
    resources:
        mem_mib = 400 * 1024,
        runtime_min = 120,
        disk_mib = 400 * 1024
    conda: '../envs/talon.yaml'
    script:
        '../scripts/talon.py'
    # shell:
    #     'talon --db {input.db} --build {params.genome_name} --threads {threads} --o {params.output_prefix} --f {input.config} > {log} 2>&1'

rule talon_abundance:
    input:
        db = 'results/talon/talon.db',
        log = 'results/talon/talon_QC.log' # ensures talon has finished
    output:
        'results/talon/abundance.tsv'
    log:
        'logs/talon/abundance.log'
    params:
        genome_name = config['genome_name'],
        annot_name = config['annot_name']
    threads: 32
    conda: '../envs/talon.yaml'
    shell:
        'talon_abundance --db {input.db} -annot {params.annot_name} --build {params.genome_name} --o results/talon > {log} 2>&1'

rule talon_filter_transcripts:
    input:
        db = 'results/talon/talon.db',
        log = 'results/talon/talon_QC.log' # ensures talon has finished
    output:
        tsv = 'results/talon/filtered_talon.tsv'
    log:
        'logs/talon/filter.log'
    params:
        annot_name = config['annot_name'],
        maxFracA = 0.5, # TODO: Move to config
        minCount = 5 # TODO: Move to config
    threads: 32
    conda: '../envs/talon.yaml'
    shell:
        'talon_filter_transcripts --db {input.db} -annot {params.annot_name} -maxFracA {params.maxFracA} -minCount {params.minCount} --o {output.tsv} > {log} 2>&1'

rule talon_create_GTF:
    input:
        db = 'results/talon/talon.db',
        log = 'results/talon/talon_QC.log' # ensures talon has finished
    output:
        # TODO: one gtf per sample/tissue? use -d option with a dataset file from a new rule
        gtf = 'results/talon/talon.gtf'
    log:
        'logs/talon/create_gtf.log'
    params:
        genome_name = config['genome_name'],
        annot_name = config['annot_name'],
        output_prefix = 'results/talon/talon'
    threads: 32
    conda: '../envs/talon.yaml'
    shell:
        'talon_create_GTF --db {input.db} -annot {params.annot_name} --build {params.genome_name} --observed --o {params.output_prefix} > {log} 2>&1'

rule talon_create_adata:
    input:
        db = 'results/talon/talon.db',
        log = 'results/talon/talon_QC.log' # ensures talon has finished
    output:
        h5ad = 'results/talon/talon.h5ad'
    log:
        'logs/talon/create_adata.log'
    params:
        genome_name = config['genome_name'],
        annot_name = config['annot_name']
    threads: 32
    conda: '../envs/talon.yaml'
    shell:
        'talon_create_adata --db {input.db} -annot {params.annot_name} --build {params.genome_name} --o {output.h5ad} > {log} 2>&1'
