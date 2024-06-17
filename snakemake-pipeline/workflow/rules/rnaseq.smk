rule rnaseq:
    input:
        expand(
            "results/rnaseq/depth_{side}_{tool}_{tissue}.txt",
            side=["inside", "outside"],
            tool=["flair", "isotools"],
            tissue=util.tissues,
        ),


rule rnaseq_depth:
    input:
        bams=lambda wildcards: util.rnaseq_alignments_for_tissue(wildcards.tissue),
        bais=lambda wildcards: expand(
            "{bam}.bai", bam=util.rnaseq_alignments_for_tissue(wildcards.tissue)
        ),
        bed="results/tss_{side}_{tool}_{tissue}.bed",
    output:
        "results/rnaseq/depth_{side}_{tool}_{tissue}.txt",
    log:
        "logs/rnaseq/depth_{side}_{tool}_{tissue}.log",
    shell:
        "samtools depth -a -b {input.bed} -o {output} {input.bams} > {log} 2>&1"
