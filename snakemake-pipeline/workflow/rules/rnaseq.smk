rule rnaseq:
    input:
        expand("results/rnaseq/depth_{tissue}.txt", tissue=util.tissues)

rule rnaseq_depth:
    input:
        bams=lambda wildcards: util.rnaseq_alignments_for_tissue(wildcards.tissue)
    output:
        "results/rnaseq/depth_{tissue}.txt"
    log:
        "logs/rnaseq/depth_{tissue}.log"
    shell:
        "(samtools depth {input} > {output}) > {log} 2>&1"
