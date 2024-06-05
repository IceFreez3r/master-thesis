import pandas as pd
import os
import logging
from isotools import Transcriptome
from isotools import __version__ as isotools_version

# set up logging
logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.INFO, filename=snakemake.log[0])
logger = logging.getLogger("isotools")
logger.info(f"This is isotools version {isotools_version}")

genome_path = snakemake.params.reference_fa
annotation_gff = snakemake.params.annotation_gff

samples = snakemake.params.samples
tissues = snakemake.params.tissues
alignments = snakemake.input.bams
metadata_path = snakemake.params.sample_table
metadata = pd.read_csv(metadata_path, sep="\t")

isoseq: Transcriptome = Transcriptome.from_reference(annotation_gff, progress_bar=False)

for sample, alignment in zip(samples, alignments):
    if not os.path.exists(alignment):
        logger.error(f"File {alignment} does not exist")
        continue
    group = metadata[metadata["sample ID"] == sample]["group"].values[0]
    isoseq.add_sample_from_bam(fn=alignment, sample_name=sample, group=group, progress_bar=False)

# compute qc metrics
isoseq.add_qc_metrics(genome_path, progress_bar=False)
# add ORF predictions
isoseq.add_orf_prediction(genome_path, progress_bar=False)

# Add tissue specific filters
group_idx = {
    gn: [i for i, sample in enumerate(isoseq.samples) if sample in grp]
    for gn, grp in isoseq.groups().items()
}
for tissue in tissues:
    filter_name = "IN" + tissue.upper()
    tissue_index = group_idx[tissue]
    expression = f"g.coverage[{tissue_index},trid].sum() >= {snakemake.params.coverage_threshold}"
    isoseq.add_filter(
        tag=filter_name, expression=expression, context="transcript", update=True
    )
    logger.info(f"Added filter {filter_name} for tissue {tissue}: {expression}")

isoseq.add_filter(
    tag="HIGH_COVER",
    expression="g.coverage.sum(0)[trid] >= 7",
    context="transcript",
    update=True,
)
isoseq.add_filter(
    tag="PERMISSIVE",
    expression="FSM or not (RTTS or INTERNAL_PRIMING or FRAGMENT)",
    context="transcript",
    update=True,
)
isoseq.add_filter(
    tag="BALANCED",
    expression="FSM or (HIGH_COVER and not (RTTS or FRAGMENT or INTERNAL_PRIMING))",
    context="transcript",
    update=True,
)
isoseq.add_filter(
    tag="STRICT",
    expression="SUBSTANTIAL and (FSM or not (RTTS or FRAGMENT or INTERNAL_PRIMING))",
    context="transcript",
    update=True,
)

isoseq.save(snakemake.output.pkl)
isoseq.write_gtf(snakemake.output.gtf, min_coverage=5, gzip=False, query="")
