import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import json
import logging
import itertools


logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.INFO, filename=snakemake.log[0])
logger = logging.getLogger("plots")


TOOLS = snakemake.params['tools']
TOOL_NAMES = snakemake.params['tool_names']
GROUPS = snakemake.params['groups']
GROUP_NAMES = snakemake.params['group_names']
CLASSIFATIONS = snakemake.params['classifications']
OUTPUT_DIR = snakemake.params['output_dir']
os.makedirs(OUTPUT_DIR, exist_ok=True)
PLOT_TITELS = snakemake.params['plot_titles']
DPI = snakemake.params['dpi']
TSS_CMAP = snakemake.params['tss_cmap']
PAS_CMAP = snakemake.params['pas_cmap']

def get_classification(group, tool) -> pd.DataFrame:
    return pd.read_csv(CLASSIFATIONS[tool][group], sep='\t')

logger.info('Importing classifications')
all = pd.DataFrame()
for (group, group_name) in zip(GROUPS, GROUP_NAMES):
    for (tool, tool_name) in zip(TOOLS, TOOL_NAMES):
        df = get_classification(group, tool)
        df.insert(0, 'group', group_name)
        df.insert(1, 'tool', tool_name)
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
# Barplot for the number of isoforms for each tool and group
agg_all = all.groupby(['tool', 'group']).agg({'TSS ratio': 'sum', 'CAGE support': 'sum', 'polyA site': 'sum', 'polyA motif': 'sum', 'start both': 'sum', 'end both': 'sum', 'isoform': 'count'}).reset_index()
agg_all['count'] = agg_all['isoform']

agg_by_category = all.groupby(['tool', 'group', 'category']).agg({'TSS ratio': 'sum', 'CAGE support': 'sum', 'polyA site': 'sum', 'polyA motif': 'sum', 'start both': 'sum', 'end both': 'sum', 'isoform': 'count'}).reset_index()
agg_by_category['count'] = agg_by_category['isoform']

agg_by_subcategory = all.groupby(['tool', 'group', 'category', 'subcategory']).agg({'TSS ratio': 'sum', 'CAGE support': 'sum', 'polyA site': 'sum', 'polyA motif': 'sum', 'start both': 'sum', 'end both': 'sum', 'isoform': 'count'}).reset_index()
agg_by_subcategory['count'] = agg_by_subcategory['isoform']

ax1 = sns.barplot(x='tool', y='count', hue='group', data=agg_all, palette='viridis')
# ax1.set_title('Number of isoforms by tool and group')
plt.savefig(os.path.join(OUTPUT_DIR, 'transcript_counts.png'), dpi=DPI)
plt.close()

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

plt.savefig(os.path.join(OUTPUT_DIR, 'transcript_counts_category.png'), dpi=DPI)
plt.close()

# ## Subcategory counts
# Counts per subcategory
df = agg_by_subcategory.groupby(['tool', 'group', 'subcategory']).agg({'count': 'sum'}).reset_index()
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

plt.savefig(os.path.join(OUTPUT_DIR, 'transcript_counts_subcategory.png'), bbox_inches='tight', dpi=DPI)
plt.close()

# ## ISM subcategory counts

df = agg_by_subcategory.loc[agg_by_subcategory['category'] == 'ISM']
# Counts per subcategory
ax1 = sns.barplot(x='tool', y='count', hue='subcategory', data=df, palette='viridis')
plt.savefig(os.path.join(OUTPUT_DIR, 'transcript_counts_subcategory_ISM.png'), dpi=DPI)
plt.close()

logger.info('Category counts plotted')

#region Heatmaps
def heatmap(df: pd.DataFrame, column, export_name, header_suffix='', **params):
    df = df.copy()
    df.loc[:,'relative_metric'] = df[column] / df['count']
    df.loc[:,'annotation'] = df['relative_metric'].map('{:,.1%}'.format) + \
                            '\n(' + df[column].astype(str) + '/' + df['count'].astype(str) + ')'
    # Reshape the data using pivot
    heatmap_data = df.pivot(index='tool', columns='group', values='relative_metric')
    # get max value for each group
    max_values = heatmap_data.max(axis=0)

    # Annotate each cell with the numeric value and the count
    annot = df.pivot(index='tool', columns='group', values='annotation')

    # Plot the heatmap
    plt.figure(figsize=(3 * heatmap_data.shape[1], 1 + heatmap_data.shape[0]))
    ax = sns.heatmap(heatmap_data, vmin=0, vmax=1, annot=annot, fmt='', linewidths=0.5, **params)
    # Make max values bold
    for i in range(heatmap_data.shape[0]):
        for j in range(heatmap_data.shape[1]):
            if heatmap_data.iloc[i, j] == max_values.iloc[j]:
                ax.texts[i * len(heatmap_data.columns) + j].set_fontweight('bold')

    if PLOT_TITELS:
        plt.title(f'Heatmap of {column}{header_suffix}')
    plt.xlabel('')
    plt.ylabel('')
    output_path = os.path.join(OUTPUT_DIR, f'{export_name}.png')
    plt.savefig(output_path, bbox_inches='tight', dpi=DPI)
    logger.info(f'Heatmap of {column}{header_suffix} plotted and saved to {output_path}')
    plt.close()

def heatmaps(df: pd.DataFrame, column, filters: list[list[tuple[str, pd.Series]]], **params):
    '''
    :param df: DataFrame to filter
    :param column: Column to plot
    :param filters: List of list of tuples with the name and a filter in form of a pandas Series
        Passed to itertools.product to generate all combinations. Use an empty string and None to skip
    :param params: Additional parameters to pass to the heatmap function
    '''
    for filter_tuple in itertools.product(*filters):
        df_filtered = df.copy(deep=True)
        suffix = ''
        for name, f in filter_tuple:
            if name == '' and f is None:
                continue
            df_filtered = df_filtered.loc[f]
            suffix += ' ' + name
        df_filtered = df_filtered.groupby(['tool', 'group']).sum().reset_index()
        heatmap(df_filtered, column, (column + suffix).replace(' ', '_'), suffix, **params)

mono_exons = agg_by_subcategory.loc[agg_by_subcategory['subcategory'] == 'mono-exon'].groupby(['group', 'tool']).sum().reset_index()
heatmap(mono_exons, 'CAGE support', header_suffix=' for Monoexons', export_name='CAGE_support_monoexons', cmap=TSS_CMAP)

FSM = agg_by_subcategory['category'] == 'FSM'
ISM = agg_by_subcategory['category'] == 'ISM'
NIC = agg_by_subcategory['category'] == 'NIC'
NNC = agg_by_subcategory['category'] == 'NNC'

no_monoexon = agg_by_subcategory['subcategory'] != 'mono-exon'
no_3prime = agg_by_subcategory['subcategory'] != '3prime_fragment'
no_5prime = agg_by_subcategory['subcategory'] != '5prime_fragment'

heatmaps(agg_by_subcategory, 'CAGE support', [
    [('', None), ('FSM', FSM), ('ISM', ISM), ('NIC', NIC), ('NNC', NNC)],
    [('', None), ('no monoexons', no_monoexon)],
    [('', None), ('no 3prime', no_3prime)]
    ], cmap=TSS_CMAP)
heatmaps(agg_by_subcategory, 'TSS ratio', [
    [('', None), ('FSM', FSM), ('ISM', ISM), ('NIC', NIC), ('NNC', NNC)],
    [('', None), ('no monoexons', no_monoexon)],
    [('', None), ('no 3prime', no_3prime)]
    ], cmap=TSS_CMAP)
heatmaps(agg_by_subcategory, 'polyA site', [
    [('', None), ('FSM', FSM), ('ISM', ISM), ('NIC', NIC), ('NNC', NNC)],
    [('', None), ('no monoexons', no_monoexon)],
    [('', None), ('no 5prime', no_5prime)]
    ], cmap=PAS_CMAP)
heatmaps(agg_by_subcategory, 'polyA motif', [
    [('', None), ('FSM', FSM), ('ISM', ISM), ('NIC', NIC), ('NNC', NNC)],
    [('', None), ('no monoexons', no_monoexon)],
    [('', None), ('no 5prime', no_5prime)]
    ], cmap=PAS_CMAP)

logger.info('Subcategory heatmaps plotted')
#endregion

#region Correlation
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
#endregion
