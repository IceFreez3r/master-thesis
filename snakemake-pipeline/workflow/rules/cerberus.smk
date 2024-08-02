import pandas as pd

localrules:
    cerberus_agg_ends_config,
    cerberus_agg_ics_config,
    # These often fail on systems without memory overcommitting
    cerberus_agg_ics,
    cerberus_agg_ends,

rule cerberus:
    input:
        annotated_gtfs = expand('results/cerberus/replace_gtf_ids/{tissue}/{tool}.gtf', tool=WORKING_TOOLS, tissue=util.tissues),


rule cerberus_reference_to_bed:
    input:
        gtf = "resources/annotation.gtf",
    output:
        bed = 'results/cerberus/triplet/{mode}/reference.bed',
    log: 'logs/cerberus/gtf_to_bed/{mode}/reference.log'
    conda:
        'cerberus'
    resources:
        mem_mb = 128 * 1024,
        runtime_min = 30,
    shell:
        'cerberus gtf_to_bed --gtf {input.gtf} --mode {wildcards.mode} -o {output.bed} > {log} 2>&1'


rule cerberus_reference_to_ics:
    input:
        gtf = "resources/annotation.gtf",
    output:
        ics = 'results/cerberus/triplet/ics/reference.ics',
    log: 'logs/cerberus/gtf_to_ics/reference.log'
    conda:
        'cerberus'
    resources:
        mem_mb = 128 * 1024,
        runtime_min = 30,
    shell:
        'cerberus gtf_to_ics --gtf {input.gtf} -o {output.ics} > {log} 2>&1'


rule cerberus_gtf_to_bed:
    input:
        gtf="results/{tool}/transcriptome/{tissue}.gtf",
    output:
        bed = 'results/cerberus/triplet/{mode}/{tissue}/{tool}.bed'
    log: 'logs/cerberus/gtf_to_bed/{mode}/{tissue}/{tool}.log'
    conda:
        'cerberus'
    resources:
        mem_mb = 128 * 1024,
        runtime_min = 30,
    shell:
        'cerberus gtf_to_bed --gtf {input.gtf} --mode {wildcards.mode} -o {output.bed} > {log} 2>&1'


rule cerberus_gtf_to_ics:
    input:
        gtf="results/{tool}/transcriptome/{tissue}.gtf",
    output:
        ics = 'results/cerberus/triplet/ics/{tissue}/{tool}.ics',
    log: 'logs/cerberus/gtf_to_ics/{tissue}/{tool}.log'
    conda:
        'cerberus'
    resources:
        mem_mb = 128 * 1024,
        runtime_min = 30,
    shell:
        'cerberus gtf_to_ics --gtf {input.gtf} -o {output.ics} > {log} 2>&1'


rule cerberus_agg_ends_config:
    '''Creates a headerless csv file with BED file path, Reference (bool), Add ends (Bool), Source name'''
    input:
        ref = 'results/cerberus/triplet/{mode}/reference.bed',
        bed = expand('results/cerberus/triplet/{{mode}}/{{tissue}}/{tool}.bed', tool=WORKING_TOOLS),
    output:
        config = 'results/cerberus/configs/agg_ends_{mode}/{tissue}.csv',
    params:
        annot_source = config["annot_source"]
    run:
        input_files = [{'BED file path': input.ref, 'Reference': True, 'Add ends': True, 'Source name': params.annot_source}]
        for bed in input.bed:
            tool = os.path.splitext(os.path.basename(bed))[0]
            input_files.append({'BED file path': bed, 'Reference': False, 'Add ends': True, 'Source name': tool})
        df = pd.DataFrame(input_files)
        df.to_csv(output.config, index=False, header=False)


rule cerberus_agg_ends:
    input:
        ref = 'results/cerberus/triplet/{mode}/reference.bed',
        bed = expand('results/cerberus/triplet/{{mode}}/{{tissue}}/{tool}.bed', tool=WORKING_TOOLS),
        config = 'results/cerberus/configs/agg_ends_{mode}/{tissue}.csv',
    output:
        bed = 'results/cerberus/agg_{mode}/{tissue}.bed'
    log: 'logs/cerberus/agg_{mode}/{tissue}.log'
    conda:
        'cerberus'
    threads: 8 # I think it uses threading under the hood
    resources:
        mem_mb = 128 * 1024,
        runtime_min = 120,
    shell:
        'cerberus agg_ends --input {input.config} --mode {wildcards.mode} -o {output.bed} > {log} 2>&1'


rule cerberus_agg_ics_config:
    '''Creates a headerless csv file with BED file path, Reference (bool), Source name'''
    input:
        ref = 'results/cerberus/triplet/ics/reference.ics',
        ics = expand('results/cerberus/triplet/ics/{{tissue}}/{tool}.ics', tool=WORKING_TOOLS),
    output:
        config = 'results/cerberus/configs/agg_ics/{tissue}.csv'
    params:
        annot_source = config["annot_source"]
    run:
        input_files = [{'BED file path': input.ref, 'Reference': True, 'Source name': params.annot_source}]
        for ics in input.ics:
            tool = os.path.splitext(os.path.basename(ics))[0]
            input_files.append({'BED file path': ics, 'Reference': False, 'Source name': tool})
        df = pd.DataFrame(input_files)
        df.to_csv(output.config, index=False, header=False)


rule cerberus_agg_ics:
    input:
        ref = 'results/cerberus/triplet/ics/reference.ics',
        ics = expand('results/cerberus/triplet/ics/{{tissue}}/{tool}.ics', tool=WORKING_TOOLS),
        config = 'results/cerberus/configs/agg_ics/{tissue}.csv',
    output:
        tsv = 'results/cerberus/agg_ics/{tissue}.tsv'
    log: 'logs/cerberus/agg_ics/{tissue}.log'
    conda:
        'cerberus'
    threads: 8 # I think it uses threading under the hood
    resources:
        mem_mb = 128 * 1024,
        runtime_min = 120,
    shell:
        'cerberus agg_ics --input {input.config} -o {output.tsv} > {log} 2>&1'


rule cerberus_write_reference:
    input:
        tss = 'results/cerberus/agg_tss/{tissue}.bed',
        tes = 'results/cerberus/agg_tes/{tissue}.bed',
        ics = 'results/cerberus/agg_ics/{tissue}.tsv',
    output:
        h5 = 'results/cerberus/reference_h5/{tissue}.h5'
    log: 'logs/cerberus/write_reference/{tissue}.log'
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
        h5 = 'results/cerberus/reference_h5/{tissue}.h5',
    output:
        h5 = 'results/cerberus/annotated/{tissue}/{tool}.h5',
    log: 'logs/cerberus/annotate/{tissue}/{tool}.log'
    conda:
        'cerberus'
    shell:
        'cerberus annotate_transcriptome --gtf {input.gtf} --h5 {input.h5} --source {wildcards.tool} -o {output.h5} > {log} 2>&1'


rule cerberus_replace_gtf_ids:
    input:
        h5 = 'results/cerberus/annotated/{tissue}/{tool}.h5',
        gtf = 'results/{tool}/transcriptome/{tissue}.gtf',
    output:
        gtf = 'results/cerberus/replace_gtf_ids/{tissue}/{tool}.gtf',
    log: 'logs/cerberus/replace_gtf_ids/{tissue}/{tool}.log'
    conda:
        'cerberus'
    resources:
        mem_mb = 128 * 1024,
        runtime_min = 30,
    shell:
        'cerberus replace_gtf_ids --h5 {input.h5} --gtf {input.gtf} --source {wildcards.tool} --update_ends --collapse -o {output.gtf} > {log} 2>&1'
