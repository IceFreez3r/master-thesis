import pandas as pd
import os
import logging
from isotools import Transcriptome
from isotools import __version__ as isotools_version

# set up logging
logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.INFO, filename=snakemake.log[0])
logger = logging.getLogger("isotools")
logger.info(f"This is isotools version {isotools_version}")

genome_path = snakemake.input.ref_fa
annotation_gff = snakemake.input.annotation_gff
metadata_path = snakemake.input.sample_table
alignments = snakemake.input.bams

group = snakemake.wildcards.group
samples = snakemake.params.samples
metadata = pd.read_csv(metadata_path, sep="\t")

extra = snakemake.params.extra
extra_add_sample_from_bam = extra.get("add_sample_from_bam", {})
extra_add_qc_metrics = extra.get("add_qc_metrics", {})
extra_write_gtf = extra.get("write_gtf")
assert "query" in extra_write_gtf, "query is required in extra.write_gtf. If you explicitly want to set it to None, use `query: ''`"

logger.info(f"Creating isotools object for group {group} from reference")
isoseq: Transcriptome = Transcriptome.from_reference(annotation_gff, progress_bar=False)

for sample, alignment in zip(samples, alignments):
    if not os.path.exists(alignment):
        logger.error(f"File {alignment} does not exist")
        continue
    group = metadata[metadata["sample ID"] == sample]["group"].values[0]
    logger.info(f"Adding sample {sample} from {alignment} to isotools object")
    isoseq.add_sample_from_bam(fn=alignment, sample_name=sample, group=group, progress_bar=False, **extra_add_sample_from_bam)

logger.info("Computing qc metrics")
isoseq.add_qc_metrics(genome_path, progress_bar=False, **extra_add_qc_metrics)

logger.info("Exporting isotools to GTF and pickle")
isoseq.write_gtf(snakemake.output.gtf, gzip=False, **extra_write_gtf)
isoseq.save(snakemake.output.pkl)
