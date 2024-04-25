rule cerberus_gtf_to_bed:
    input:
        gtf = 'data/genes.gtf'
    output:
        bed = 'results/cerberus/genes_{mode}.bed'
    conda:
        'cerberus'
    shell:
        'cerberus gtf_to_bed --gtf {input.gtf} --mode {wildcards.mode} -o {output.bed}'

rule cerberus_gtf_to_ics:
    input:
        gtf = 'data/genes.gtf'
    output:
        ics = 'results/cerberus/genes.ics'
    conda:
        'cerberus'
    shell:
        'cerberus gtf_to_ics --gtf {input.gtf} -o {output.ics}'

rule cerberus_agg_ends_config:
    '''Creates a headerless csv file with BED file path, Reference (bool), Add ends (Bool), Source name'''
    input:
        bed = 'results/cerberus/genes_{mode}.bed'
    output:
        config = 'results/cerberus/agg_ends_config.csv'
    run:
        ""

rule cerberus_agg_ends:
    input:
        config = 'results/cerberus/agg_ends_config.csv',
    output:
        bed = 'results/cerberus/agg_ends_{mode}.bed'
    conda:
        'cerberus'
    shell:
        'cerberus agg_ends --input {input.config} --mode {wildcards.mode} -o {output.bed}'

rule cerberus_agg_ics_config:
    '''Creates a headerless csv file with BED file path, Reference (bool), Add ends (Bool), Source name'''
    input:
        bed = 'results/cerberus/agg_ends_{mode}.bed'
    output:
        config = 'results/cerberus/agg_ends_config.csv'
    run:
        ""

rule cerberus_agg_ics:
    input:
        csv = 'results/cerberus/agg_ends_config.csv',
    output:
        tsv = 'results/cerberus/agg_ends_ics.tsv'
    conda:
        'cerberus'
    shell:
        'cerberus agg_ics --input {input.csv} -o {output.tsv}'

rule cerberus_write_reference:
    input:
        tss = 'results/cerberus/agg_ends_tss.bed',
        tes = 'results/cerberus/agg_ends_tes.bed',
        ics = 'results/cerberus/agg_ends_ics.tsv'
    output:
        h5 = 'results/cerberus/reference.h5'
    conda:
        'cerberus'
    shell:
        'cerberus write_reference --tss {input.tss} --tes {input.tes} --ics {input.ics} -o {output.h5}'

rule cerberus_annotate_transcriptome:
    input:
        gtf = 'data/genes.gtf',
        reference = 'results/cerberus/reference.h5'
    output:
        bed = 'results/cerberus/genes_annotated.bed' # TODO
    conda:
        'cerberus'
    shell:
        'cerberus annotate_transcriptome --gtf {input.gtf} --h5 {input.reference} -o {output.bed}'

rule cerberus_replace_ab_ids:
    input:
        reference = 'results/cerberus/reference.h5',
        talon_abundance = 'results/talon/abundance.tsv',
        source = "?"
    output:
        tsv = 'results/cerberus/abundance.tsv'
    conda:
        'cerberus'
    shell:
        'cerberus replace_ab_ids --h5 {input.reference} --ab {input.talon_abundance} --source {input.source} --collapse -o {output.tsv}'

rule cerberus_replace_gtf_ids:
    input:
        reference = 'results/cerberus/reference.h5',
        talon_gtf = 'results/talon/genes.gtf',
        source = "?"
    output:
        gtf = 'results/cerberus/genes.gtf'
    conda:
        'cerberus'
    shell:
        'cerberus replace_gtf_ids --h5 {input.reference} --gtf {input.talon_gtf} --source {input.source} --update_ends --collapse -o {output.gtf}'
