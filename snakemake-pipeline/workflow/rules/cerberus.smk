import pandas as pd

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
        bed = expand('results/cerberus/triplet/{{mode}}/{sample}.bed', sample=samples)
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
        bed = expand('results/cerberus/triplet/{{mode}}/{sample}.bed', sample=samples),
        config = 'results/cerberus/configs/agg_ends_{mode}.csv'
    output:
        bed = 'results/cerberus/agg_{mode}.bed'
    log: 'logs/cerberus/agg_ends/{mode}.log'
    conda:
        'cerberus'
    shell:
        'cerberus agg_ends --input {input.config} --mode {wildcards.mode} -o {output.bed} > {log} 2>&1'

rule cerberus_agg_ics_config:
    '''Creates a headerless csv file with BED file path, Reference (bool), Add ends (Bool), Source name'''
    input:
        ref = 'results/cerberus/reference_ics.ics',
        ics = expand('results/cerberus/triplet/ics/{sample}.ics', sample=samples)
    output:
        config = 'results/cerberus/configs/agg_ics.csv'
    run:
        df = pd.DataFrame({'BED file path': [input.ref], 'Reference': [True], 'Source name': ['ENCODE']})
        for ics in input.ics:
            df = pd.concat([df, pd.DataFrame({'BED file path': [ics], 'Reference': [False], 'Source name': ['TALON']})])
        df.to_csv(output.config, index=False, header=False)

rule cerberus_agg_ics:
    input:
        ref = 'results/cerberus/reference_ics.ics',
        ics = expand('results/cerberus/triplet/ics/{sample}.ics', sample=samples),
        config = 'results/cerberus/configs/agg_ics.csv'
    output:
        tsv = 'results/cerberus/agg_ics.tsv'
    log: 'logs/cerberus/agg_ics.log'
    conda:
        'cerberus'
    shell:
        'cerberus agg_ics --input {input.config} -o {output.tsv} > {log} 2>&1'

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
    shell:
        'cerberus write_reference --tss {input.tss} --tes {input.tes} --ics {input.ics} -o {output.h5} > {log} 2>&1'

rule cerberus_annotate_transcriptome:
    input:
        gtf = 'data/genes.gtf',
        reference = 'results/cerberus/reference.h5'
    output:
        bed = 'results/cerberus/genes_annotated.bed' # TODO
    log: 'logs/cerberus/annotate_transcriptome.log'
    conda:
        'cerberus'
    shell:
        'cerberus annotate_transcriptome --gtf {input.gtf} --h5 {input.reference} -o {output.bed} > {log} 2>&1'

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

rule cerberus_replace_gtf_ids:
    input:
        reference = 'results/cerberus/reference.h5',
        talon_gtf = 'results/talon/genes.gtf',
    output:
        gtf = 'results/cerberus/genes.gtf'
    params:
        source = "TALON"
    log: 'logs/cerberus/replace_gtf_ids.log'
    conda:
        'cerberus'
    shell:
        'cerberus replace_gtf_ids --h5 {input.reference} --gtf {input.talon_gtf} --source {params.source} --update_ends --collapse -o {output.gtf} > {log} 2>&1'
