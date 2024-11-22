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

tissue = snakemake.wildcards.tissue
samples = snakemake.params.samples
query = snakemake.params.query
metadata = pd.read_csv(metadata_path, sep="\t")

unify_ends = snakemake.params.get("unify_ends")

logger.info(f"Creating isotools object for tissue {tissue} from reference")
isoseq: Transcriptome = Transcriptome.from_reference(annotation_gff, progress_bar=False)

for sample, alignment in zip(samples, alignments):
    if not os.path.exists(alignment):
        logger.error(f"File {alignment} does not exist")
        continue
    group = metadata[metadata["sample ID"] == sample]["group"].values[0]
    logger.info(f"Adding sample {sample} from {alignment} to isotools object")
    isoseq.add_sample_from_bam(fn=alignment, sample_name=sample, group=group, progress_bar=False)

logger.info("Computing qc metrics. Unifying ends: {unify_ends}")
isoseq.add_qc_metrics(genome_path, progress_bar=False, unify_ends=unify_ends)

logger.info("Exporting isotools to GTF and pickle")
isoseq.write_gtf(snakemake.output.gtf, gzip=False, query=query)
isoseq.save(snakemake.output.pkl)
