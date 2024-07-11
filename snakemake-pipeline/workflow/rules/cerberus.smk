import pandas as pd

localrules:
    cerberus_agg_ends_config,
    cerberus_agg_ics_config,
    gz_gtf,
    # These often fail on systems without memory overcommitting
    cerberus_agg_ics,
    cerberus_agg_ends,

rule cerberus:
    input:
        annotated_gtfs = expand('results/cerberus/{tool}/replace_gtf_ids/{tissue}.gtf', tool=WORKING_TOOLS, tissue=util.tissues),

# rule filter_gtf:
#     '''Filters out `source` attribute from the gtf. Would throw an error in cerberus_gtf_to_bed and is redundant anyways.'''
#     input:
#         gtf = 'resources/transcriptome/{tissue}.gtf'
#     output:
#         gtf = 'resources/transcriptome/{tissue}_filtered.gtf'
#     run:
#         with open(input.gtf, 'r') as f, open(output.gtf, 'w') as out:
#             for line in f:
#                 attr = line.split('\t')[-1].split(';')
#                 attr = [a for a in attr if 'source' not in a]
#                 out.write('\t'.join(line.split('\t')[:-1] + [';'.join(attr)]))

rule cerberus_reference_to_bed:
    input:
        gtf = "resources/annotation.gtf"
    output:
        bed = 'results/cerberus/reference_{mode}.bed'
    log: 'logs/cerberus/reference/to_bed_{mode}.log'
    conda:
        'cerberus'
    resources:
        mem_mb = 128 * 1024,
        runtime_min = 120,
    shell:
        'cerberus gtf_to_bed --gtf {input.gtf} --mode {wildcards.mode} -o {output.bed} > {log} 2>&1'

rule cerberus_reference_to_ics:
    input:
        gtf = "resources/annotation.gtf"
    output:
        ics = 'results/cerberus/reference_ics.ics'
    log: 'logs/cerberus/reference/to_ics.log'
    conda:
        'cerberus'
    resources:
        mem_mb = 128 * 1024,
        runtime_min = 120,
    shell:
        'cerberus gtf_to_ics --gtf {input.gtf} -o {output.ics} > {log} 2>&1'

rule cerberus_gtf_to_bed:
    input:
        gtf="results/{tool}/transcriptome/{tissue}.gtf",
    output:
        bed = 'results/cerberus/{tool}/triplet/{mode}/{tissue}.bed'
    log: 'logs/cerberus/{tool}/gtf_to_bed/{mode}_{tissue}.log'
    conda:
        'cerberus'
    resources:
        mem_mb = 128 * 1024,
        runtime_min = 120,
    shell:
        'cerberus gtf_to_bed --gtf {input.gtf} --mode {wildcards.mode} -o {output.bed} > {log} 2>&1'

rule cerberus_gtf_to_ics:
    input:
        gtf="results/{tool}/transcriptome/{tissue}.gtf",
    output:
        ics = 'results/cerberus/{tool}/triplet/ics/{tissue}.ics'
    log: 'logs/cerberus/{tool}/gtf_to_ics/{tissue}.log'
    conda:
        'cerberus'
    resources:
        mem_mb = 128 * 1024,
        runtime_min = 120,
    shell:
        'cerberus gtf_to_ics --gtf {input.gtf} -o {output.ics} > {log} 2>&1'

rule cerberus_agg_ends_config:
    '''Creates a headerless csv file with BED file path, Reference (bool), Add ends (Bool), Source name'''
    input:
        ref = 'results/cerberus/reference_{mode}.bed',
        bed = expand('results/cerberus/{{tool}}/triplet/{{mode}}/{tissue}.bed', tissue=util.tissues)
    output:
        config = 'results/cerberus/{tool}/configs/agg_ends_{mode}.csv'
    run:
        df = pd.DataFrame({'BED file path': [input.ref], 'Reference': [True], 'Add ends': [True], 'Source name': ['ENCODE']})
        for bed in input.bed:
            df = pd.concat([df, pd.DataFrame({'BED file path': [bed], 'Reference': [False], 'Add ends': [True], 'Source name': [wildcards.tool]})])
        df.to_csv(output.config, index=False, header=False)

rule cerberus_agg_ends:
    input:
        ref = 'results/cerberus/reference_{mode}.bed',
        bed = expand('results/cerberus/{{tool}}/triplet/{{mode}}/{tissue}.bed', tissue=util.tissues),
        config = 'results/cerberus/{tool}/configs/agg_ends_{mode}.csv'
    output:
        bed = 'results/cerberus/{tool}/agg_{mode}.bed'
    log: 'logs/cerberus/{tool}/agg_ends/{mode}.log'
    conda:
        'cerberus'
    resources:
        mem_mb = 128 * 1024,
        runtime_min = 120,
    shell:
        'cerberus agg_ends --input {input.config} --mode {wildcards.mode} -o {output.bed} > {log} 2>&1'

rule cerberus_agg_ics_config:
    '''Creates a headerless csv file with BED file path, Reference (bool), Add ends (Bool), Source name'''
    input:
        ref = 'results/cerberus/reference_ics.ics',
        ics = expand('results/cerberus/{{tool}}/triplet/ics/{tissue}.ics', tissue=util.tissues),
    output:
        config = 'results/cerberus/{tool}/configs/agg_ics.csv'
    run:
        df = pd.DataFrame({'BED file path': [input.ref], 'Reference': [True], 'Source name': ['ENCODE']})
        for ics in input.ics:
            df = pd.concat([df, pd.DataFrame({'BED file path': [ics], 'Reference': [False], 'Source name': ['TALON']})])
        df.to_csv(output.config, index=False, header=False)

rule cerberus_agg_ics:
    input:
        ref = 'results/cerberus/reference_ics.ics',
        ics = expand('results/cerberus/{{tool}}/triplet/ics/{tissue}.ics', tissue=util.tissues),
        config = 'results/cerberus/{tool}/configs/agg_ics.csv'
    output:
        tsv = 'results/cerberus/{tool}/agg_ics.tsv'
    log: 'logs/cerberus/{tool}/agg_ics.log'
    conda:
        'cerberus'
    threads: 8 # I think it uses threading under the hood
    resources:
        mem_mb = 128 * 1024,
        runtime_min = 120,
    shell:
        'cerberus agg_ics --input {input.config} -o {output.tsv} > {log} 2>&1'

# rule cerberus_agg_ics_script:
#     input:
#         ref = 'results/cerberus/reference_ics.ics',
#         ics = expand('results/cerberus/{{tool}}/triplet/ics/{tissue}.ics', tissue=util.tissues),
#     output:
#         tsv = 'results/cerberus/{tool}/agg_ics.tsv',
#         config = 'results/cerberus/{tool}/configs/agg_ics.csv'
#     params:
#         limit = config["cerberus"]['aggregate_at_once_limit']
#     log: 'logs/cerberus/{tool}/agg_ics.log'
#     conda:
#         'cerberus'
#     benchmark:
#         'benchmarks/cerberus/{tool}/agg_ics.tsv'
#     threads: 8 # I think it uses threading under the hood
#     resources:
#         mem_mb = 128 * 1024,
#         runtime_min = 120,
#     script:
#         '../scripts/cerberus_agg_ics.py'

rule cerberus_write_reference:
    input:
        tss = 'results/cerberus/{tool}/agg_tss.bed',
        tes = 'results/cerberus/{tool}/agg_tes.bed',
        ics = 'results/cerberus/{tool}/agg_ics.tsv'
    output:
        h5 = 'results/cerberus/{tool}/reference.h5'
    log: 'logs/cerberus/{tool}/write_reference.log'
    conda:
        'cerberus'
    threads: 16 # I think it uses threading under the hood
    resources:
        mem_mb = 128 * 1024,
        runtime_min = 120,
    shell:
        'cerberus write_reference --tss {input.tss} --tes {input.tes} --ics {input.ics} -o {output.h5} > {log} 2>&1'

rule cerberus_annotate_transcriptome:
    input:
        gtf = "results/{tool}/transcriptome/{tissue}.gtf",
        h5 = 'results/cerberus/{tool}/reference.h5'
    output:
        gtf = 'results/cerberus/{tool}/annotated/{tissue}.gtf'
    log: 'logs/cerberus/{tool}/annotate/{tissue}.log'
    params:
        source = "TALON"
    conda:
        'cerberus'
    shell:
        'cerberus annotate_transcriptome --gtf {input.gtf} --h5 {input.h5} --source {params.source} -o {output.gtf} > {log} 2>&1'

rule cerberus_replace_ab_ids:
    input:
        reference = 'results/cerberus/{tool}/reference.h5',
        talon_abundance = 'results/talon/abundance.tsv',
    output:
        tsv = 'results/cerberus/{tool}/abundance.tsv'
    params:
        source = "TALON"
    log: 'logs/cerberus/{tool}/replace_ab_ids.log'
    conda:
        'cerberus'
    resources:
        mem_mb = 128 * 1024,
        runtime_min = 120,
    shell:
        'cerberus replace_ab_ids --h5 {input.reference} --ab {input.talon_abundance} --source {params.source} --collapse -o {output.tsv} > {log} 2>&1'

rule gz_gtf:
    input:
        gtf = 'results/cerberus/{tool}/annotated/{tissue}.gtf'
    output:
        gtf_gz = 'results/cerberus/{tool}/annotated/{tissue}.gtf.gz'
    shell:
        'gzip -c {input.gtf} > {output.gtf_gz}'

# TODO: Unclear if it wants a .gz or not
rule cerberus_replace_gtf_ids:
    input:
        reference = 'results/cerberus/{tool}/reference.h5',
        gtf_gz = 'results/cerberus/{tool}/annotated/{tissue}.gtf.gz',
    output:
        gtf = 'results/cerberus/{tool}/replace_gtf_ids/{tissue}.gtf'
    params:
        source = "TALON"
    log: 'logs/cerberus/{tool}/replace_gtf_ids/{tissue}.log'
    conda:
        'cerberus'
    resources:
        mem_mb = 128 * 1024,
        runtime_min = 120,
    shell:
        'cerberus replace_gtf_ids --h5 {input.reference} --gtf {input.gtf_gz} --source {params.source} --update_ends --collapse -o {output.gtf} > {log} 2>&1'
