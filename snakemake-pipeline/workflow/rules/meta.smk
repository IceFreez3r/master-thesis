rule dag:
    output:
        "results/meta/dag.pdf"
    shell:
        "snakemake --forceall --rulegraph | dot -Tpdf > {output}"
