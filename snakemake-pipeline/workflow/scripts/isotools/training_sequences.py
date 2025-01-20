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
pkl_path = snakemake.input.pkl
sqanti_annotation = snakemake.input.sqanti_annotation

output_prefix = snakemake.params.output_prefix
extra = snakemake.params.extra
logger.info(extra)

logger.info(f"Importing isotools object from pickle {pkl_path}")
isoseq: Transcriptome = Transcriptome.load(pkl_path)
logger.info(f"Importing SQANTI annotation from {sqanti_annotation}")
isoseq.import_sqanti_classification(sqanti_annotation, progress_bar=False)
logger.info(f"Exporting end sequences to {output_prefix}")
isoseq.export_end_sequences(genome_path, output_prefix, **extra)
logger.info("Done")
