import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import json
import logging


logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.INFO, filename=snakemake.log[0])
logger = logging.getLogger("plots")


TISSUES = snakemake.params['tissues']
TOOLS = snakemake.params['tools']
TOOLNAMES = snakemake.params['toolnames']
print(snakemake.params['classifications_dict'])
CLASSIFATIONS = json.load(snakemake.params['classifications_dict'])
OUTPUT_DIR = snakemake.params['output_dir']
PLOT_TITELS = snakemake.params['plot_titles']

def get_classification(tissue, tool) -> pd.DataFrame:
    return pd.read_csv(CLASSIFATIONS[tissue][tool], sep='\t')

logger.info('Importing classifications')
all = pd.DataFrame()
for tissue in TISSUES:
    for (tool, name) in zip(TOOLS, TOOLNAMES):
        df = get_classification(tissue, tool)
        df.insert(0, 'tissue', tissue)
        df.insert(1, 'tool', name)
        # Rename and make boolean
        df['category'] = df['structural_category']
        df['TSS ratio'] = df['ratio_TSS'] > 1.5
        df['CAGE support'] = df['within_CAGE_peak']
        df['polyA site'] = df['within_polyA_site']
        df['polyA motif'] = df['polyA_motif_found']
        df['start both'] = df['TSS ratio'] & df['CAGE support']
        df['end both'] = df['polyA site'] & df['polyA motif']
        all = pd.concat([all, df], ignore_index=True)

# Rename categories
all['category'] = all['category'].replace({
    'intergenic': 'Intergenic',
    'antisense': 'Antisense',
    'genic': 'Genic',
    'genic_intron': 'Genic Intron',
    'fusion': 'Fusion',
    'full-splice_match': 'FSM',
    'incomplete-splice_match': 'ISM',
    'novel_in_catalog': 'NIC',
    'novel_not_in_catalog': 'NNC'
})
logger.info('Classifications imported')

# # Transcript Counts
# Barplot for the number of isoforms for each tool and tissue
agg_all = all.groupby(['tool', 'tissue']).agg({'TSS ratio': 'sum', 'CAGE support': 'sum', 'polyA site': 'sum', 'polyA motif': 'sum', 'start both': 'sum', 'end both': 'sum', 'isoform': 'count'}).reset_index()
agg_all['count'] = agg_all['isoform']

agg_by_category = all.groupby(['tool', 'tissue', 'category']).agg({'TSS ratio': 'sum', 'CAGE support': 'sum', 'polyA site': 'sum', 'polyA motif': 'sum', 'start both': 'sum', 'end both': 'sum', 'isoform': 'count'}).reset_index()
agg_by_category['count'] = agg_by_category['isoform']

agg_by_subcategory = all.groupby(['tool', 'tissue', 'category', 'subcategory']).agg({'TSS ratio': 'sum', 'CAGE support': 'sum', 'polyA site': 'sum', 'polyA motif': 'sum', 'start both': 'sum', 'end both': 'sum', 'isoform': 'count'}).reset_index()
agg_by_subcategory['count'] = agg_by_subcategory['isoform']

ax1 = sns.barplot(x='tool', y='count', hue='tissue', data=agg_all, palette='viridis')
# ax1.set_title('Number of isoforms by tool and tissue')
plt.savefig(os.path.join(OUTPUT_DIR, 'transcript_counts.png'))

# # Category Counts
# Counts per category
# ax1 = sns.barplot(x='tool', y='count', hue='category', data=agg_by_category, palette='viridis')
# no_stringtie = ~agg_by_category['tool'].isin(['stringtie'])
# ax1 = sns.barplot(x='tool', y='count', hue='category', data=agg_by_category[no_stringtie], palette='viridis')

# Based on https://matplotlib.org/stable/gallery/subplots_axes_and_figures/broken_axis.html
f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': (1, 3)})
sns.barplot(x='tool', y='count', hue='category', data=agg_by_category, palette='viridis', ax=ax1, legend=False)
# ax1.set_title('Number of isoforms by tool and category')
ax1.set_ylabel('')
sns.barplot(x='tool', y='count', hue='category', data=agg_by_category, palette='viridis', ax=ax2)
ax2.legend(loc='upper center', bbox_to_anchor=(0.3, 1.45), ncol=2)

ax1.set_ylim(245_000, 265_000)  # outliers only
ax2.set_ylim(0, 55_000)  # most of the data

# hide the spines between ax and ax2
ax1.spines.bottom.set_visible(False)
ax2.spines.top.set_visible(False)
ax1.xaxis.tick_top()
ax1.tick_params(labeltop=False)  # don't put tick labels at the top
ax2.xaxis.tick_bottom()

d = .5  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)

plt.savefig(os.path.join(OUTPUT_DIR, 'transcript_counts_category.png'))

# ## Subcategory counts
# Counts per subcategory
# df = agg_by_subcategory.groupby(['tool', 'tissue', 'subcategory']).agg({'count': 'sum'}).reset_index()
# ax1 = sns.barplot(x='tool', y='count', hue='subcategory', data=df, palette='viridis')
# no_stringtie = ~df['tool'].isin(['stringtie'])
# ax1 = sns.barplot(x='tool', y='count', hue='subcategory', data=df[no_stringtie], palette='viridis')

# split y between 40k and 140k
# move legend left outside the plot
f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, height_ratios=(1, 5), figsize=(12, 6))
sns.barplot(x='tool', y='count', hue='subcategory', data=df, palette='viridis', ax=ax1, legend=False)
sns.barplot(x='tool', y='count', hue='subcategory', data=df, palette='viridis', ax=ax2)

ax1.set_ylabel('')

ax1.set_ylim(151_000, 163_000)  # outliers only
ax2.set_ylim(0, 63_000)  # most of the data

# hide the spines between ax and ax2
ax1.spines.bottom.set_visible(False)
ax2.spines.top.set_visible(False)
ax1.xaxis.tick_top()
ax1.tick_params(labeltop=False)  # don't put tick labels at the top
ax2.xaxis.tick_bottom()

d = .5  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)

# Move legend outside the plot
ax2.legend(loc='upper left', bbox_to_anchor=(1, 1.3), ncol=1)

plt.savefig(os.path.join(OUTPUT_DIR, 'transcript_counts_subcategory.png'), bbox_inches='tight')

# ## ISM subcategory counts

df = agg_by_subcategory.loc[agg_by_subcategory['category'] == 'ISM']
# Counts per subcategory
ax1 = sns.barplot(x='tool', y='count', hue='subcategory', data=df, palette='viridis')
plt.savefig(os.path.join(OUTPUT_DIR, 'transcript_counts_subcategory_ISM.png'))

logger.info('Category counts plotted')
# # All

def heatmap(df: pd.DataFrame, column, export_name, header_suffix='', **params):
    df = df.copy()
    df.loc[:,'relative_metric'] = df[column] / df['count']
    df.loc[:,'annotation'] = df['relative_metric'].map('{:,.1%}'.format) + \
                            '\n(' + df[column].astype(str) + '/' + df['count'].astype(str) + ')'
    # Reshape the data using pivot
    heatmap_data = df.pivot(index='tool', columns='tissue', values='relative_metric')
    # get max value for each tissue
    max_values = heatmap_data.max(axis=0)

    # Annotate each cell with the numeric value and the count
    annot = df.pivot(index='tool', columns='tissue', values='annotation')

    # Plot the heatmap
    plt.figure(figsize=(12, 1 + len(df['tool'].unique())))
    ax = sns.heatmap(heatmap_data, vmin=0, vmax=1, annot=annot, fmt='', linewidths=0.5, **params)
    # Make max values bold
    for i in range(heatmap_data.shape[0]):
        for j in range(heatmap_data.shape[1]):
            if heatmap_data.iloc[i, j] == max_values.iloc[j]:
                ax.texts[i * len(heatmap_data.columns) + j].set_fontweight('bold')

    if PLOT_TITELS:
        plt.title(f'Heatmap of {column}{header_suffix}')
    plt.xlabel('Tissue')
    plt.ylabel('Tool')
    plt.savefig(os.path.join(OUTPUT_DIR, f'{export_name}.png'), bbox_inches='tight')

heatmap(agg_all, 'CAGE support', export_name='CAGE_support', cmap='YlGnBu')
heatmap(agg_all, 'TSS ratio', export_name='TSS_ratio', cmap='YlGnBu')
# heatmap(agg_all, 'start both', cmap='YlGnBu')

heatmap(agg_all, 'polyA site', export_name='PolyA_site', cmap='magma_r')
heatmap(agg_all, 'polyA motif', export_name='PolyA_motif', cmap='magma_r')
# heatmap(agg_all, 'end both', cmap='magma_r')

# # By Category

non_fsm_df = agg_by_category.loc[agg_by_category['category'] != 'FSM'].groupby(['tissue', 'tool']).sum().reset_index()

# ## Starts

heatmap(agg_by_category.loc[agg_by_category['category'] == 'FSM'], 'CAGE support', ' for FSM', export_name='CAGE_support_FSM', cmap='YlGnBu')
heatmap(agg_by_category.loc[agg_by_category['category'] == 'ISM'], 'CAGE support', ' for ISM', export_name='CAGE_support_ISM', cmap='YlGnBu')
heatmap(agg_by_category.loc[agg_by_category['category'] == 'NIC'], 'CAGE support', ' for NIC', export_name='CAGE_support_NIC', cmap='YlGnBu')
heatmap(agg_by_category.loc[agg_by_category['category'] == 'NNC'], 'CAGE support', ' for NNC', export_name='CAGE_support_NNC', cmap='YlGnBu')
heatmap(non_fsm_df, 'CAGE support', ' for non-FSM', export_name='CAGE_support_non_FSM', cmap='YlGnBu')

# heatmap(agg_by_category.loc[agg_by_category['category'] == 'FSM'], 'TSS ratio', ' for FSM', cmap='YlGnBu')
# heatmap(agg_by_category.loc[agg_by_category['category'] == 'ISM'], 'TSS ratio', ' for ISM', cmap='YlGnBu')
# heatmap(non_fsm_df, 'TSS ratio', ' for non-FSM', cmap='YlGnBu')

# heatmap(agg_by_category.loc[agg_by_category['category'] == 'FSM'], 'start both', ' for FSM', cmap='YlGnBu')
# heatmap(agg_by_category.loc[agg_by_category['category'] == 'ISM'], 'start both', ' for ISM', cmap='YlGnBu')
# heatmap(non_fsm_df, 'start both', ' for non-FSM', cmap='YlGnBu')

# ## Ends

heatmap(agg_by_category.loc[agg_by_category['category'] == 'FSM'], 'polyA site', ' for FSM', export_name="polyA_site_FSM", cmap='magma_r')
heatmap(agg_by_category.loc[agg_by_category['category'] == 'ISM'], 'polyA site', ' for ISM', export_name="polyA_site_ISM", cmap='magma_r')
heatmap(agg_by_category.loc[agg_by_category['category'] == 'NIC'], 'polyA site', ' for NIC', export_name="polyA_site_NIC", cmap='magma_r')
heatmap(agg_by_category.loc[agg_by_category['category'] == 'NNC'], 'polyA site', ' for NNC', export_name="polyA_site_NNC", cmap='magma_r')
# heatmap(non_fsm_df, 'polyA site', ' for non-FSM', cmap='magma_r')

heatmap(agg_by_category.loc[agg_by_category['category'] == 'FSM'], 'polyA motif', ' for FSM', export_name="polyA_motif_FSM", cmap='magma_r')
heatmap(agg_by_category.loc[agg_by_category['category'] == 'ISM'], 'polyA motif', ' for ISM', export_name="polyA_motif_ISM", cmap='magma_r')
heatmap(agg_by_category.loc[agg_by_category['category'] == 'NIC'], 'polyA motif', ' for NIC', export_name="polyA_motif_NIC", cmap='magma_r')
heatmap(agg_by_category.loc[agg_by_category['category'] == 'NNC'], 'polyA motif', ' for NNC', export_name="polyA_motif_NNC", cmap='magma_r')
# heatmap(non_fsm_df, 'polyA motif', ' for non-FSM', cmap='magma_r')

# heatmap(agg_by_category.loc[agg_by_category['category'] == 'FSM'], 'end both', ' for FSM', cmap='magma_r')
# heatmap(agg_by_category.loc[agg_by_category['category'] == 'ISM'], 'end both', ' for ISM', cmap='magma_r')
# heatmap(non_fsm_df, 'end both', ' for non-FSM', cmap='magma_r')

logger.info('Category heatmaps plotted')

# # Subcategory

# ## Only Monoexons

mono_exons = agg_by_subcategory.loc[agg_by_subcategory['subcategory'] == 'mono-exon'].groupby(['tissue', 'tool']).sum().reset_index()
heatmap(mono_exons, 'CAGE support', ' for Monoexons', export_name='CAGE_support_monoexons', cmap='YlGnBu')
# heatmap(mono_exons, 'TSS ratio', ' for Monoexons', cmap='YlGnBu')
# heatmap(mono_exons, 'start both', ' for Monoexons', cmap='YlGnBu')

# ## ISM w/o Monoexons

no_mono_ISM = agg_by_subcategory.loc[(agg_by_subcategory['subcategory'] != 'mono-exon') & (agg_by_subcategory['category'] == 'ISM')].groupby(['tissue', 'tool']).sum().reset_index()
heatmap(no_mono_ISM, 'CAGE support', ' for ISM w/o Monoexons', export_name='CAGE_support_ISM_no_monoexons', cmap='YlGnBu')
# heatmap(no_mono_ISM, 'TSS ratio', ' for ISM w/o Monoexons', cmap='YlGnBu')
# heatmap(no_mono_ISM, 'start both', ' for ISM w/o Monoexons', cmap='YlGnBu')

# # All w/o Monoexons

no_mono_all = agg_by_subcategory.loc[(agg_by_subcategory['subcategory'] != 'mono-exon')].groupby(['tissue', 'tool']).sum().reset_index()
heatmap(no_mono_all, 'CAGE support', ' w/o Monoexons', export_name='CAGE_support_no_monoexons', cmap='YlGnBu')
# heatmap(no_mono_all, 'TSS ratio', ' w/o Monoexons', cmap='YlGnBu')
# heatmap(no_mono_all, 'start both', ' w/o Monoexons', cmap='YlGnBu')

# # No monoexons and 3' fragments

heatmap(agg_by_subcategory.loc[(agg_by_subcategory['subcategory'] != 'mono-exon') & (agg_by_subcategory['subcategory'] != '3prime_fragment')].groupby(['tissue', 'tool']).sum().reset_index(),
        'CAGE support',
        ' w/o Monoexons and 3\' Fragments',
        export_name='CAGE_support_no_monoexons_no_3fragment',
        cmap='YlGnBu')
heatmap(agg_by_subcategory.loc[(agg_by_subcategory['subcategory'] != 'mono-exon') & (agg_by_subcategory['subcategory'] != '3prime_fragment') & (agg_by_subcategory['category'] == 'FSM')].groupby(['tissue', 'tool']).sum().reset_index(),
        'CAGE support',
        ' w/o Monoexons and 3\' Fragments',
        export_name='CAGE_support_FSM_no_monoexons_no_3fragment',
        cmap='YlGnBu')
heatmap(agg_by_subcategory.loc[(agg_by_subcategory['subcategory'] != 'mono-exon') & (agg_by_subcategory['subcategory'] != '3prime_fragment') & (agg_by_subcategory['category'] == 'ISM')].groupby(['tissue', 'tool']).sum().reset_index(),
        'CAGE support',
        ' w/o Monoexons and 3\' Fragments',
        export_name='CAGE_support_ISM_no_monoexons_no_3fragment',
        cmap='YlGnBu')
heatmap(agg_by_subcategory.loc[(agg_by_subcategory['subcategory'] != 'mono-exon') & (agg_by_subcategory['subcategory'] != '3prime_fragment') & (agg_by_subcategory['category'] == 'NIC')].groupby(['tissue', 'tool']).sum().reset_index(),
        'CAGE support',
        ' w/o Monoexons and 3\' Fragments',
        export_name='CAGE_support_NIC_no_monoexons_no_3fragment',
        cmap='YlGnBu')
heatmap(agg_by_subcategory.loc[(agg_by_subcategory['subcategory'] != 'mono-exon') & (agg_by_subcategory['subcategory'] != '3prime_fragment') & (agg_by_subcategory['category'] == 'NNC')].groupby(['tissue', 'tool']).sum().reset_index(),
        'CAGE support',
        ' w/o Monoexons and 3\' Fragments',
        export_name='CAGE_support_NNC_no_monoexons_no_3fragment',
        cmap='YlGnBu')

logger.info('Subcategory heatmaps plotted')

# # Correlation
STATS_FILE = snakemake.output['stats']

with open(STATS_FILE, 'w') as f:
    def correlation(df: pd.DataFrame, column1, column2):
        df = df.copy()
        df = df.loc[df['count'] > 0]
        # Calculate the correlation
        corr, p = pearsonr(df[column1], df[column2])
        f.write(f'Correlation between {column1} and {column2}: {corr:.3f} (p={p:.2e})\n')

    correlation(agg_all, 'CAGE support', 'TSS ratio')
    correlation(agg_all, 'polyA motif', 'polyA site')

logger.info('Correlations calculated')
logger.info('Done')

