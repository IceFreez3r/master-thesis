import pyranges as pr
import sys


with open(snakemake.log[0], "w") as log:
    sys.stderr = sys.stdout = log
    # Filter out all lines from the gtf without strand info
    tr = pr.read_gtf(snakemake.input.gtf)
    tr = tr[(tr.Strand == "+") | (tr.Strand == "-")]
    tr.to_gtf(snakemake.output.gtf)
