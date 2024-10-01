import pyranges as pr
import sys


with open(snakemake.log[0], "w") as log:
    sys.stderr = sys.stdout = log
    # Filter out all lines from the gtf without strand info
    tr = pr.read_gtf(snakemake.input.gtf)
    num_tr = len(tr)
    tr = tr[(tr.Strand == "+") | (tr.Strand == "-")]
    num_tr_filtered = len(tr)
    log.write(f"Filtered out {num_tr - num_tr_filtered} transcripts without strand info\n")
    tr.to_gtf(snakemake.output.gtf)
