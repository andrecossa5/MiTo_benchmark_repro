"""
scLT systems tree comparison.
"""

import os
from mito_utils.utils import *
from mito_utils.plotting_base import *
matplotlib.use('macOSX')


##


# Get metrics
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro'
path_data = os.path.join(path_main, 'results', 'MI_TO_bench', 'scLT_comparison')
path_results = os.path.join(path_main, 'results', 'MI_TO_bench', 'scLT_comparison')

# Read tables
df_metrics = pd.read_csv(os.path.join(path_data, 'scLT_comparison_tree_metrics.csv'), index_col=0)

# Select metrics to plot
metrics = [
    'n_cells', 'n_characters', 'median_char_per_cell', 'density', 'median_CI', 'median_RI',
    'corr_distances', 'median_support', 'median_support_biggest_clades',
    'frac_remaining_nodes_after_mutationless_edges_collapse'
]
df = df_metrics[['method']+metrics]
names = dict(zip(
    metrics, 
    ['n cells', 'n characters', 'n var alleles per cell', 'Density',
     'CI', 'RI', 'Tree vs char- dists correlation',
     'Support', 'Support biggest clades', '% nodes after tree pruning'
    ]
))


##


# Dataset choice
fig, ax = plt.subplots(figsize=(4,4))
df_ = df['method'].value_counts().to_frame('n samples')
bar(df_, 'n samples', ax=ax, a=.7, s=.75, edgecolor='k')
format_ax(xlabel='scLT system', ylabel='n samples', xticks=df_.index, ax=ax, reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'scLT_system_n_samples.png'), dpi=700)


##


# Viz
fig, axs = plt.subplots(2,5,figsize=(10,5),sharex=True)

order = df.groupby('method')['median_support_biggest_clades'].median().sort_values().index
colors = create_palette(df, 'method', 'Spectral_r')

for i,metric in enumerate(metrics):
    ax = axs.ravel()[i]
    box(df, x='method', y=metric, ax=ax, c=colors, order=order)
    format_ax(ax=ax, reduced_spines=True, ylabel=names[metric], rotx=90)

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'scLT_system_overall_comparison.png'), dpi=700)


##


