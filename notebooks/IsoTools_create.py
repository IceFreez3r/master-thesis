# %%
# %load_ext autoreload
# %autoreload 2

# %%
import sys
sys.path.append('/home/lankenau/isotools/src')

# %%
import pandas as pd
import matplotlib.pyplot as plt
import os
import math
import copy

path = '/project/hfa_work/ENCODE/data'
alignment_path = 'alignment_v45'
genome_file = 'GRCh38.p14.genome.fa'
genome_path = os.path.join(path, 'gencode_human/version_45', genome_file)

# %%
metadata_file = 'reads/metadata_tissue.tsv'
metadata = pd.read_csv(os.path.join(path, metadata_file), sep='\t')
metadata

# %%
import logging
from isotools import Transcriptome
from isotools import __version__ as isotools_version
# set up logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
logger = logging.getLogger('isotools')
logger.info(f'This is isotools version {isotools_version}')

# %%
tissues = metadata['group'].unique()
tissues

# %%
annotation_file = os.path.join(path, 'gencode_human/version_45', 'gencode.v45.chr_patch_hapl_scaff.annotation_sorted.gff3.gz')
# create one IsoTools transcriptome object from the reference annotation per tissue
isoseq: Transcriptome = Transcriptome.from_reference(annotation_file)

# %%
for i, row in metadata.iterrows():
    sample_name = row['sample ID']
    # file is the full (wrong) path, we just need the filename without the extension
    sample_file = os.path.join(path, alignment_path, row['file'].split('/')[-1].split('.')[0] + '_aligned.bam')
    if not os.path.exists(sample_file):
        logger.error(f'File {sample_file} does not exist')
        continue
    group = row['group']
    isoseq.add_sample_from_bam(fn=sample_file, sample_name=sample_name, group=group)
isoseq.sample_table

# %%
# compute qc metrics
isoseq.add_qc_metrics(genome_path)
# add ORF predictions
isoseq.add_orf_prediction(genome_path)

# %%
group_idx = {gn:[i for i,sa in enumerate(isoseq.samples) if sa in grp] for gn,grp in isoseq.groups().items()}
for tissue in tissues:
    filter_name = 'IN' + tissue.upper()
    tissue_index = group_idx[tissue]
    # TODO: Lower coverage threshold, high for now to reduce time
    expression = f'g.coverage[{tissue_index},trid].sum() >= 10'
    isoseq.add_filter(tag=filter_name, expression=expression, context='transcript', update=True)
    print(f'Added filter {filter_name} for tissue {tissue}: {expression}')

# %%
isoseq.add_filter(tag='HIGH_COVER',
                  expression='g.coverage.sum(0)[trid] >= 7',
                  context='transcript',
                  update=True)
isoseq.add_filter(tag='PERMISSIVE',
                  expression='FSM or not (RTTS or INTERNAL_PRIMING or FRAGMENT)',
                  context='transcript',
                  update=True)
isoseq.add_filter(tag='BALANCED',
                  expression='FSM or (HIGH_COVER and not (RTTS or FRAGMENT or INTERNAL_PRIMING))',
                  context='transcript',
                  update=True)
isoseq.add_filter(tag='STRICT',
                  expression='SUBSTANTIAL and (FSM or not (RTTS or FRAGMENT or INTERNAL_PRIMING))',
                  context='transcript',
                  update=True)

# %%
isoseq.save('results/isoseq_v45.pkl')

# %%
isoseq.write_gtf('results/isoseq_v45.gtf', min_coverage=5, gzip=False, query="")

# %%
isoseq.groups()


