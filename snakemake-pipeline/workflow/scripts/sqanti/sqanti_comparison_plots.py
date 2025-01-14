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
# # All

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
    plt.savefig(os.path.join(OUTPUT_DIR, f'{export_name}.png'), bbox_inches='tight', dpi=DPI)
    plt.close()

heatmap(agg_all, 'CAGE support', export_name='CAGE_support', cmap=TSS_CMAP)
heatmap(agg_all, 'TSS ratio', export_name='TSS_ratio', cmap=TSS_CMAP)
# heatmap(agg_all, 'start both', cmap=TSS_CMAP)

heatmap(agg_all, 'polyA site', export_name='PolyA_site', cmap=PAS_CMAP)
heatmap(agg_all, 'polyA motif', export_name='PolyA_motif', cmap=PAS_CMAP)
# heatmap(agg_all, 'end both', cmap=PAS_CMAP)

# # By Category

non_fsm_df = agg_by_category.loc[agg_by_category['category'] != 'FSM'].groupby(['group', 'tool']).sum().reset_index()

# ## Starts

heatmap(agg_by_category.loc[agg_by_category['category'] == 'FSM'], 'CAGE support',
        header_suffix=' for FSM', export_name='CAGE_support_FSM', cmap=TSS_CMAP)
heatmap(agg_by_category.loc[agg_by_category['category'] == 'ISM'], 'CAGE support',
        header_suffix=' for ISM', export_name='CAGE_support_ISM', cmap=TSS_CMAP)
heatmap(agg_by_category.loc[agg_by_category['category'] == 'NIC'], 'CAGE support',
        header_suffix=' for NIC', export_name='CAGE_support_NIC', cmap=TSS_CMAP)
heatmap(agg_by_category.loc[agg_by_category['category'] == 'NNC'], 'CAGE support',
        header_suffix=' for NNC', export_name='CAGE_support_NNC', cmap=TSS_CMAP)
heatmap(non_fsm_df, 'CAGE support',
        header_suffix=' for non-FSM', export_name='CAGE_support_non_FSM', cmap=TSS_CMAP)

# heatmap(agg_by_category.loc[agg_by_category['category'] == 'FSM'], 'TSS ratio',
#       header_suffix=' for FSM', cmap=TSS_CMAP)
# heatmap(agg_by_category.loc[agg_by_category['category'] == 'ISM'], 'TSS ratio',
#       header_suffix=' for ISM', cmap=TSS_CMAP)
# heatmap(non_fsm_df, 'TSS ratio',
#       header_suffix=' for non-FSM', cmap=TSS_CMAP)

# heatmap(agg_by_category.loc[agg_by_category['category'] == 'FSM'], 'start both',
#       header_suffix=' for FSM', cmap=TSS_CMAP)
# heatmap(agg_by_category.loc[agg_by_category['category'] == 'ISM'], 'start both',
#       header_suffix=' for ISM', cmap=TSS_CMAP)
# heatmap(non_fsm_df, 'start both',
#       header_suffix=' for non-FSM', cmap=TSS_CMAP)

# ## Ends

heatmap(agg_by_category.loc[agg_by_category['category'] == 'FSM'], 'polyA site',
        header_suffix=' for FSM', export_name="PolyA_site_FSM", cmap=PAS_CMAP)
heatmap(agg_by_category.loc[agg_by_category['category'] == 'ISM'], 'polyA site',
        header_suffix=' for ISM', export_name="PolyA_site_ISM", cmap=PAS_CMAP)
heatmap(agg_by_category.loc[agg_by_category['category'] == 'NIC'], 'polyA site',
        header_suffix=' for NIC', export_name="PolyA_site_NIC", cmap=PAS_CMAP)
heatmap(agg_by_category.loc[agg_by_category['category'] == 'NNC'], 'polyA site',
        header_suffix=' for NNC', export_name="PolyA_site_NNC", cmap=PAS_CMAP)
# heatmap(non_fsm_df, 'polyA site',
#       header_suffix=' for non-FSM', cmap=PAS_CMAP)

heatmap(agg_by_category.loc[agg_by_category['category'] == 'FSM'], 'polyA motif',
        header_suffix=' for FSM', export_name="PolyA_motif_FSM", cmap=PAS_CMAP)
heatmap(agg_by_category.loc[agg_by_category['category'] == 'ISM'], 'polyA motif',
        header_suffix=' for ISM', export_name="PolyA_motif_ISM", cmap=PAS_CMAP)
heatmap(agg_by_category.loc[agg_by_category['category'] == 'NIC'], 'polyA motif',
        header_suffix=' for NIC', export_name="PolyA_motif_NIC", cmap=PAS_CMAP)
heatmap(agg_by_category.loc[agg_by_category['category'] == 'NNC'], 'polyA motif',
        header_suffix=' for NNC', export_name="PolyA_motif_NNC", cmap=PAS_CMAP)
# heatmap(non_fsm_df, 'polyA motif',
#       header_suffix=' for non-FSM', cmap=PAS_CMAP)

# heatmap(agg_by_category.loc[agg_by_category['category'] == 'FSM'], 'end both',
#       header_suffix=' for FSM', cmap=PAS_CMAP)
# heatmap(agg_by_category.loc[agg_by_category['category'] == 'ISM'], 'end both',
#       header_suffix=' for ISM', cmap=PAS_CMAP)
# heatmap(non_fsm_df, 'end both', header_suffix=' for non-FSM', cmap=PAS_CMAP)

logger.info('Category heatmaps plotted')

# # Subcategory

# ## Only Monoexons

mono_exons = agg_by_subcategory.loc[agg_by_subcategory['subcategory'] == 'mono-exon'].groupby(['group', 'tool']).sum().reset_index()
heatmap(mono_exons, 'CAGE support', header_suffix=' for Monoexons',
        export_name='CAGE_support_monoexons', cmap=TSS_CMAP)
# heatmap(mono_exons, 'TSS ratio', header_suffix=' for Monoexons', cmap=TSS_CMAP)
# heatmap(mono_exons, 'start both', header_suffix=' for Monoexons', cmap=TSS_CMAP)

# ## ISM w/o Monoexons

no_mono_ISM = agg_by_subcategory.loc[(agg_by_subcategory['subcategory'] != 'mono-exon') & (agg_by_subcategory['category'] == 'ISM')].groupby(['group', 'tool']).sum().reset_index()
heatmap(no_mono_ISM, 'CAGE support', header_suffix=' for ISM w/o Monoexons',
        export_name='CAGE_support_ISM_no_monoexons', cmap=TSS_CMAP)
# heatmap(no_mono_ISM, 'TSS ratio', header_suffix=' for ISM w/o Monoexons', cmap=TSS_CMAP)
# heatmap(no_mono_ISM, 'start both', header_suffix=' for ISM w/o Monoexons', cmap=TSS_CMAP)

# # All w/o Monoexons

no_mono_all = agg_by_subcategory.loc[(agg_by_subcategory['subcategory'] != 'mono-exon')].groupby(['group', 'tool']).sum().reset_index()
heatmap(no_mono_all, 'CAGE support', header_suffix=' w/o Monoexons',
        export_name='CAGE_support_no_monoexons', cmap=TSS_CMAP)
# heatmap(no_mono_all, 'TSS ratio', header_suffix=' w/o Monoexons', cmap=TSS_CMAP)
# heatmap(no_mono_all, 'start both', header_suffix=' w/o Monoexons', cmap=TSS_CMAP)

# # No monoexons and 3' fragments

heatmap(agg_by_subcategory.loc[(agg_by_subcategory['subcategory'] != 'mono-exon') & (agg_by_subcategory['subcategory'] != '3prime_fragment')].groupby(['group', 'tool']).sum().reset_index(),
        'CAGE support',
        header_suffix=' w/o Monoexons and 3\' Fragments',
        export_name='CAGE_support_no_monoexons_no_3fragment',
        cmap=TSS_CMAP)
heatmap(agg_by_subcategory.loc[(agg_by_subcategory['subcategory'] != 'mono-exon') & (agg_by_subcategory['subcategory'] != '3prime_fragment') & (agg_by_subcategory['category'] == 'FSM')].groupby(['group', 'tool']).sum().reset_index(),
        'CAGE support',
        header_suffix=' w/o Monoexons and 3\' Fragments',
        export_name='CAGE_support_FSM_no_monoexons_no_3fragment',
        cmap=TSS_CMAP)
heatmap(agg_by_subcategory.loc[(agg_by_subcategory['subcategory'] != 'mono-exon') & (agg_by_subcategory['subcategory'] != '3prime_fragment') & (agg_by_subcategory['category'] == 'ISM')].groupby(['group', 'tool']).sum().reset_index(),
        'CAGE support',
        header_suffix=' w/o Monoexons and 3\' Fragments',
        export_name='CAGE_support_ISM_no_monoexons_no_3fragment',
        cmap=TSS_CMAP)
heatmap(agg_by_subcategory.loc[(agg_by_subcategory['subcategory'] != 'mono-exon') & (agg_by_subcategory['subcategory'] != '3prime_fragment') & (agg_by_subcategory['category'] == 'NIC')].groupby(['group', 'tool']).sum().reset_index(),
        'CAGE support',
        header_suffix=' w/o Monoexons and 3\' Fragments',
        export_name='CAGE_support_NIC_no_monoexons_no_3fragment',
        cmap=TSS_CMAP)
heatmap(agg_by_subcategory.loc[(agg_by_subcategory['subcategory'] != 'mono-exon') & (agg_by_subcategory['subcategory'] != '3prime_fragment') & (agg_by_subcategory['category'] == 'NNC')].groupby(['group', 'tool']).sum().reset_index(),
        'CAGE support',
        header_suffix=' w/o Monoexons and 3\' Fragments',
        export_name='CAGE_support_NNC_no_monoexons_no_3fragment',
        cmap=TSS_CMAP)

heatmap(agg_by_subcategory.loc[(agg_by_subcategory['subcategory'] != 'mono-exon') & (agg_by_subcategory['subcategory'] != '5prime_fragment')].groupby(['group', 'tool']).sum().reset_index(),
        'polyA site',
        header_suffix=' w/o Monoexons and 5\' Fragments',
        export_name='PolyA_site_no_monoexons_no_5fragment',
        cmap=PAS_CMAP)
heatmap(agg_by_subcategory.loc[(agg_by_subcategory['subcategory'] != 'mono-exon') & (agg_by_subcategory['subcategory'] != '5prime_fragment') & (agg_by_subcategory['category'] == 'FSM')].groupby(['group', 'tool']).sum().reset_index(),
        'polyA site',
        header_suffix=' w/o Monoexons and 5\' Fragments',
        export_name='PolyA_site_FSM_no_monoexons_no_5fragment',
        cmap=PAS_CMAP)
heatmap(agg_by_subcategory.loc[(agg_by_subcategory['subcategory'] != 'mono-exon') & (agg_by_subcategory['subcategory'] != '5prime_fragment') & (agg_by_subcategory['category'] == 'ISM')].groupby(['group', 'tool']).sum().reset_index(),
        'polyA site',
        header_suffix=' w/o Monoexons and 5\' Fragments',
        export_name='PolyA_site_ISM_no_monoexons_no_5fragment',
        cmap=PAS_CMAP)
heatmap(agg_by_subcategory.loc[(agg_by_subcategory['subcategory'] != 'mono-exon') & (agg_by_subcategory['subcategory'] != '5prime_fragment') & (agg_by_subcategory['category'] == 'NIC')].groupby(['group', 'tool']).sum().reset_index(),
        'polyA site',
        header_suffix=' w/o Monoexons and 5\' Fragments',
        export_name='PolyA_site_NIC_no_monoexons_no_5fragment',
        cmap=PAS_CMAP)
heatmap(agg_by_subcategory.loc[(agg_by_subcategory['subcategory'] != 'mono-exon') & (agg_by_subcategory['subcategory'] != '5prime_fragment') & (agg_by_subcategory['category'] == 'NNC')].groupby(['group', 'tool']).sum().reset_index(),
        'polyA site',
        header_suffix=' w/o Monoexons and 5\' Fragments',
        export_name='PolyA_site_NNC_no_monoexons_no_5fragment',
        cmap=PAS_CMAP)

heatmap(agg_by_subcategory.loc[(agg_by_subcategory['subcategory'] != 'mono-exon') & (agg_by_subcategory['subcategory'] != '5prime_fragment')].groupby(['group', 'tool']).sum().reset_index(),
        'polyA motif',
        header_suffix=' w/o Monoexons and 5\' Fragments',
        export_name='PolyA_motif_no_monoexons_no_5fragment',
        cmap=PAS_CMAP)
heatmap(agg_by_subcategory.loc[(agg_by_subcategory['subcategory'] != 'mono-exon') & (agg_by_subcategory['subcategory'] != '5prime_fragment') & (agg_by_subcategory['category'] == 'FSM')].groupby(['group', 'tool']).sum().reset_index(),
        'polyA motif',
        header_suffix=' w/o Monoexons and 5\' Fragments',
        export_name='PolyA_motif_FSM_no_monoexons_no_5fragment',
        cmap=PAS_CMAP)
heatmap(agg_by_subcategory.loc[(agg_by_subcategory['subcategory'] != 'mono-exon') & (agg_by_subcategory['subcategory'] != '5prime_fragment') & (agg_by_subcategory['category'] == 'ISM')].groupby(['group', 'tool']).sum().reset_index(),
        'polyA motif',
        header_suffix=' w/o Monoexons and 5\' Fragments',
        export_name='PolyA_motif_ISM_no_monoexons_no_5fragment',
        cmap=PAS_CMAP)
heatmap(agg_by_subcategory.loc[(agg_by_subcategory['subcategory'] != 'mono-exon') & (agg_by_subcategory['subcategory'] != '5prime_fragment') & (agg_by_subcategory['category'] == 'NIC')].groupby(['group', 'tool']).sum().reset_index(),
        'polyA motif',
        header_suffix=' w/o Monoexons and 5\' Fragments',
        export_name='PolyA_motif_NIC_no_monoexons_no_5fragment',
        cmap=PAS_CMAP)
heatmap(agg_by_subcategory.loc[(agg_by_subcategory['subcategory'] != 'mono-exon') & (agg_by_subcategory['subcategory'] != '5prime_fragment') & (agg_by_subcategory['category'] == 'NNC')].groupby(['group', 'tool']).sum().reset_index(),
        'polyA motif',
        header_suffix=' w/o Monoexons and 5\' Fragments',
        export_name='PolyA_motif_NNC_no_monoexons_no_5fragment',
        cmap=PAS_CMAP)

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

