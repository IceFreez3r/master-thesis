rule star:
    input:
        expand("results/star/{sample}/{sample}.SJ.out.tab", sample=util.rnaseq_samples),
        expand("results/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam", sample=util.rnaseq_samples),

rule star_index:
    input:
        ref_fa = "resources/reference.fa",
    output:
        index = directory("resources/star_index"),
    log:
        "logs/star/index.log"
    params:
        limitGenomeGenerateRAM = lambda wc, resources: resources.mem_mb * 1024 * 1024,
        extra = config["star"]["index"]["extra"],
    threads: 64
    resources:
        mem_mb=128 * 1024,
        runtime_min=3 * 60,
        disk_mb=128 * 1024
    conda:
        "../envs/star.yaml"
    shell:
        """
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output.index} \
             --genomeFastaFiles {input.ref_fa} \
             --limitGenomeGenerateRAM {params.limitGenomeGenerateRAM} {params.extra} > {log} 2>&1
        """


rule star_map:
    input:
        index = "resources/star_index",
        fastq = lambda wildcards: util.rnaseq_reads_for_sample(wildcards.sample),
    output:
        sj_tab = "results/star/{sample}/{sample}.SJ.out.tab",
        bam = "results/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
    log:
        "results/star/{sample}/{sample}.Log.out",
        "results/star/{sample}/{sample}.Log.progress.out",
        "results/star/{sample}/{sample}.Log.final.out",
        main = "logs/star/{sample}.log",
    params:
        extra = config["star"]["align"]["extra"],
        limitBAMsortRAM = lambda wc, resources: f"--limitBAMsortRAM {resources.mem_mb * 1024 * 1024}",
        output_prefix = lambda wc, output: output.sj_tab.replace("SJ.out.tab", ""),
    threads: 16
    resources:
        mem_mb=128 * 1024,
        runtime_min=2 * 60,
        disk_mb=128 * 1024
    conda:
        "../envs/star.yaml"
    shell:
        """
        (
            STAR --runThreadN {threads} --genomeDir {input.index} --readFilesIn {input.fastq} \
                 --outFileNamePrefix {params.output_prefix} --readFilesCommand zcat \
                 --outSAMtype BAM SortedByCoordinate {params.limitBAMsortRAM} {params.extra} \
                 --genomeLoad NoSharedMemory \
        ) > {log.main} 2>&1
        """
