import pandas as pd

# localrules: cerberus_agg_ends_config, cerberus_agg_ics_config

rule cerberus:
    input:
        annotated_gtfs = expand('results/cerberus/replace_gtf_ids/{sample}.gtf', sample=util.samples),

rule filter_gtf:
    '''Filters out `source` attribute from the gtf. Would throw an error in cerberus_gtf_to_bed and is redundant anyways.'''
    input:
        gtf = 'resources/transcriptome/{sample}.gtf'
    output:
        gtf = 'resources/transcriptome/{sample}_filtered.gtf'
    run:
        with open(input.gtf, 'r') as f, open(output.gtf, 'w') as out:
            for line in f:
                attr = line.split('\t')[-1].split(';')
                attr = [a for a in attr if 'source' not in a]
                out.write('\t'.join(line.split('\t')[:-1] + [';'.join(attr)]))

rule cerberus_reference_to_bed:
    input:
        gtf = "resources/reference_annotation.gtf"
    output:
        bed = 'results/cerberus/reference_{mode}.bed'
    log: 'logs/cerberus/reference/to_bed_{mode}.log'
    conda:
        'cerberus'
    shell:
        'cerberus gtf_to_bed --gtf {input.gtf} --mode {wildcards.mode} -o {output.bed} > {log} 2>&1'

rule cerberus_reference_to_ics:
    input:
        gtf = "resources/reference_annotation.gtf"
    output:
        ics = 'results/cerberus/reference_ics.ics'
    log: 'logs/cerberus/reference/to_ics.log'
    conda:
        'cerberus'
    shell:
        'cerberus gtf_to_ics --gtf {input.gtf} -o {output.ics} > {log} 2>&1'

rule cerberus_gtf_to_bed:
    input:
        gtf = 'resources/transcriptome/{sample}_filtered.gtf'
    output:
        bed = 'results/cerberus/triplet/{mode}/{sample}.bed'
    log: 'logs/cerberus/gtf_to_bed/{mode}_{sample}.log'
    conda:
        'cerberus'
    shell:
        'cerberus gtf_to_bed --gtf {input.gtf} --mode {wildcards.mode} -o {output.bed} > {log} 2>&1'

rule cerberus_gtf_to_ics:
    input:
        gtf = 'resources/transcriptome/{sample}_filtered.gtf'
    output:
        ics = 'results/cerberus/triplet/ics/{sample}.ics'
    log: 'logs/cerberus/gtf_to_ics/{sample}.log'
    conda:
        'cerberus'
    shell:
        'cerberus gtf_to_ics --gtf {input.gtf} -o {output.ics} > {log} 2>&1'

rule cerberus_agg_ends_config:
    '''Creates a headerless csv file with BED file path, Reference (bool), Add ends (Bool), Source name'''
    input:
        ref = 'results/cerberus/reference_{mode}.bed',
        bed = expand('results/cerberus/triplet/{{mode}}/{sample}.bed', sample=util.samples)
    output:
        config = 'results/cerberus/configs/agg_ends_{mode}.csv'
    run:
        df = pd.DataFrame({'BED file path': [input.ref], 'Reference': [True], 'Add ends': [True], 'Source name': ['ENCODE']})
        for bed in input.bed:
            df = pd.concat([df, pd.DataFrame({'BED file path': [bed], 'Reference': [False], 'Add ends': [True], 'Source name': ['TALON']})])
        df.to_csv(output.config, index=False, header=False)

rule cerberus_agg_ends:
    input:
        ref = 'results/cerberus/reference_{mode}.bed',
        bed = expand('results/cerberus/triplet/{{mode}}/{sample}.bed', sample=util.samples),
        config = 'results/cerberus/configs/agg_ends_{mode}.csv'
    output:
        bed = 'results/cerberus/agg_{mode}.bed'
    log: 'logs/cerberus/agg_ends/{mode}.log'
    conda:
        'cerberus'
    shell:
        'cerberus agg_ends --input {input.config} --mode {wildcards.mode} -o {output.bed} > {log} 2>&1'

# rule cerberus_agg_ics_config:
#     '''Creates a headerless csv file with BED file path, Reference (bool), Add ends (Bool), Source name'''
#     input:
#         ref = 'results/cerberus/reference_ics.ics',
#         ics = expand('results/cerberus/triplet/ics/{sample}.ics', sample=util.samples)
#     output:
#         config = 'results/cerberus/configs/agg_ics.csv'
#     run:
#         df = pd.DataFrame({'BED file path': [input.ref], 'Reference': [True], 'Source name': ['ENCODE']})
#         for ics in input.ics:
#             df = pd.concat([df, pd.DataFrame({'BED file path': [ics], 'Reference': [False], 'Source name': ['TALON']})])
#         df.to_csv(output.config, index=False, header=False)

# rule cerberus_agg_ics:
#     input:
#         ref = 'results/cerberus/reference_ics.ics',
#         ics = expand('results/cerberus/triplet/ics/{sample}.ics', sample=util.samples),
#         config = 'results/cerberus/configs/agg_ics.csv'
#     output:
#         tsv = 'results/cerberus/agg_ics.tsv'
#     log: 'logs/cerberus/agg_ics.log'
#     conda:
#         'cerberus'
#     threads: 8 # I think it uses threading under the hood
#     resources:
#         mem_mib = 100 * 1024,
#         runtime_min = 120,
#     shell:
#         'cerberus agg_ics --input {input.config} -o {output.tsv} > {log} 2>&1'

rule cerberus_agg_ics_script:
    input:
        ref = 'results/cerberus/reference_ics.ics',
        ics = expand('results/cerberus/triplet/ics/{sample}.ics', sample=util.samples),
    output:
        tsv = 'results/cerberus/agg_ics.tsv',
        config = 'results/cerberus/configs/agg_ics.csv'
    params:
        limit = config['aggregate_at_once_limit']
    log: 'logs/cerberus/agg_ics.log'
    conda:
        'cerberus'
    benchmark:
        'benchmarks/cerberus/agg_ics.tsv'
    threads: 8 # I think it uses threading under the hood
    resources:
        mem_mib = 900 * 1024,
        runtime_min = 120,
    script:
        '../scripts/cerberus_agg_ics.py'

rule cerberus_write_reference:
    input:
        tss = 'results/cerberus/agg_tss.bed',
        tes = 'results/cerberus/agg_tes.bed',
        ics = 'results/cerberus/agg_ics.tsv'
    output:
        h5 = 'results/cerberus/reference.h5'
    log: 'logs/cerberus/write_reference.log'
    conda:
        'cerberus'
    threads: 16 # I think it uses threading under the hood
    resources:
        mem_mib = 100 * 1024,
        runtime_min = 120,
    shell:
        'cerberus write_reference --tss {input.tss} --tes {input.tes} --ics {input.ics} -o {output.h5} > {log} 2>&1'

rule cerberus_annotate_transcriptome:
    input:
        gtf = 'resources/transcriptome/{sample}_filtered.gtf',
        h5 = 'results/cerberus/reference.h5'
    output:
        gtf = 'results/cerberus/annotated/{sample}.gtf'
    log: 'logs/cerberus/annotate/{sample}.log'
    params:
        source = "TALON"
    conda:
        'cerberus'
    shell:
        'cerberus annotate_transcriptome --gtf {input.gtf} --h5 {input.h5} --source {params.source} -o {output.gtf} > {log} 2>&1'

rule cerberus_replace_ab_ids:
    input:
        reference = 'results/cerberus/reference.h5',
        talon_abundance = 'results/talon/abundance.tsv',
    output:
        tsv = 'results/cerberus/abundance.tsv'
    params:
        source = "TALON"
    log: 'logs/cerberus/replace_ab_ids.log'
    conda:
        'cerberus'
    shell:
        'cerberus replace_ab_ids --h5 {input.reference} --ab {input.talon_abundance} --source {params.source} --collapse -o {output.tsv} > {log} 2>&1'

rule gz_gtf:
    input:
        gtf = 'results/cerberus/annotated/{sample}.gtf'
    output:
        gtf_gz = 'results/cerberus/annotated/{sample}.gtf.gz'
    shell:
        'gzip -c {input.gtf} > {output.gtf_gz}'

# TODO: Unclear if it wants a .gz or not
rule cerberus_replace_gtf_ids:
    input:
        reference = 'results/cerberus/reference.h5',
        gtf_gz = 'results/cerberus/annotated/{sample}.gtf.gz',
    output:
        gtf = 'results/cerberus/replace_gtf_ids/{sample}.gtf'
    params:
        source = "TALON"
    log: 'logs/cerberus/replace_gtf_ids/{sample}.log'
    conda:
        'cerberus'
    shell:
        'cerberus replace_gtf_ids --h5 {input.reference} --gtf {input.gtf_gz} --source {params.source} --update_ends --collapse -o {output.gtf} > {log} 2>&1'
