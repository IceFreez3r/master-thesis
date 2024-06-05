import logging
from isotools import Transcriptome
from isotools import __version__ as isotools_version

# set up logging
logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.INFO, filename=snakemake.log[0])
logger = logging.getLogger("isotools")
logger.info(f"This is isotools version {isotools_version}")

isoseq: Transcriptome = Transcriptome.load(snakemake.input.pkl)

for tissue, tissue_gtf in zip(snakemake.params.tissues, snakemake.output.tissue_gtfs):
    isoseq.write_gtf(
        f"results/isotools_{tissue}.gtf",
        min_coverage=0,
        gzip=False,
        query=f"IN{tissue.upper()} and ({snakemake.params.query})",
    )
