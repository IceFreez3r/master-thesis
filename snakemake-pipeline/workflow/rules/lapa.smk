import pandas as pd
import os

rule lapa:
    input:
        'results/lapa/links/links.tsv'

rule chromsizes:
    output:
        chromsizes = 'resources/genome.fasta.chromsizes',
    params:
        fasta = config['reference_fa']
    log: 'logs/lapa/chromsizes.log'
    shell:
        '(faidx {params.fasta} -i chromsizes > {output.chromsizes}) > {log} 2>&1'

rule gencode_utr_fix:
    input:
        'results/annotation.gtf',
    output:
        'results/annotation_utr_fix.gtf'
    log:
        'logs/lapa/gencode_utr_fix.log'
    conda:
        'lapa'
    resources:
        mem_mib = 32 * 1024,
        runtime_min = 10,
    shell:
        'gencode_utr_fix --input_gtf {input} --output_gtf {output} > {log} 2>&1'

rule lapa_config:
    input:
        bams = expand(os.path.join(config['alignment_dir'], '{sample}_aligned.bam'), sample=util.experiments),
    output:
        samples = 'results/lapa/samples.csv'
    threads: 4
    resources:
        mem_mib = 4 * 1024,
        runtime_min = 5,
    log: 'logs/lapa/config.log'
    run:
        with open(output['samples'], 'w') as f:
            df = pd.DataFrame()
            for sample in util.samples:
                df = pd.concat([df, pd.DataFrame({'sample': sample, 'dataset': util.tissue_for_sample(sample), 'path': util.alignment_for_sample(sample)}, index=[0])])
            df.to_csv(f, index=False)

rule lapa_tes:
    input:
        samples = 'results/lapa/samples.csv',
        chromsizes = 'resources/genome.fasta.chromsizes',
        gtf = 'results/annotation_utr_fix.gtf' if config['annot_source'] == "GENCODE" else 'results/annotation.gtf',
    output:
        directory('results/lapa/tes'),
        'results/lapa/tes/polyA_clusters.bed'
    params:
        fasta = config["reference_fa"],
        output_dir = 'results/lapa/tes/',
        min_replication_rate = 0.75
    log:
        'logs/lapa/tes.log'
    threads: 16
    resources:
        mem_mib = 32 * 1024,
        runtime_min = 120,
    conda:
        'lapa'
    shell:
        'lapa --alignment {input.samples} --fasta {params.fasta} --annotation {input.gtf} --chrom_sizes {input.chromsizes} --output_dir {params.output_dir} --min_replication_rate {params.min_replication_rate} > {log} 2>&1'

rule lapa_tss:
    input:
        samples = 'results/lapa/samples.csv',
        chromsizes = 'resources/genome.fasta.chromsizes',
        gtf = 'results/annotation_utr_fix.gtf' if config['annot_source'] == "GENCODE" else 'results/annotation.gtf',
    output:
        directory('results/lapa/tss'),
        'results/lapa/tss/tss_clusters.bed'
    params:
        fasta = config["reference_fa"],
        output_dir = 'results/lapa/tss/',
        min_replication_rate = 0.75
    log:
        'logs/lapa/tss.log'
    threads: 16
    resources:
        mem_mib = 32 * 1024,
        runtime_min = 120,
    conda:
        'lapa'
    shell:
        'lapa_tss --alignment {input.samples} --fasta {params.fasta} --annotation {input.gtf} --chrom_sizes {input.chromsizes} --output_dir {params.output_dir} --min_replication_rate {params.min_replication_rate} > {log} 2>&1'

rule lapa_link:
    input:
        samples = 'results/lapa/samples.csv',
        lapa_tes = 'results/lapa/tes/polyA_clusters.bed',
        lapa_tss = 'results/lapa/tss/tss_clusters.bed',
    output:
        'results/lapa/links/links.tsv'
    params:
        lapa_dir = 'results/lapa/',
        lapa_tss_dir = 'results/lapa/tss/',
    log:
        'logs/lapa/links.log'
    threads: 16
    resources:
        mem_mib = 32 * 1024,
        runtime_min = 120,
    conda:
        'lapa'
    shell:
        'lapa_link_tss_to_te --alignments {input.samples} --lapa_dir {params.lapa_dir} --lapa_tss_dir {params.lapa_tss_dir} --output {output} > {log} 2>&1'

# TODO: needed?
rule lapa_correct_talon:
    input:
    # TODO: add the correct input files
        links = 'results/lapa/links/links.tsv',
        read_annot = 'results/talon/read_annot.tsv',
        talon_gtf = 'results/talon/talon.gtf',
        talon_abundance = 'results/talon/talon_abundance.tsv',
    output:
        gtf = 'results/lapa/talon_corrected.gtf',
        abundance = 'results/lapa/talon_abundance_corrected.tsv'
    log:
        'logs/lapa/talon_corrected.log'
    conda:
        'lapa'
    shell:
        'lapa_correct_talon --links {input.links} --read_annot {input.read_annot} --gtf_input {input.talon_gtf} --abundance_input {input.talon_abundance} --gtf_output {output.gtf} --abundance_output {output.abundance} > {log} 2>&1'
