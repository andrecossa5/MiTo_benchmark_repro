"""
Supp Fig ... preprocessing pipelines comparison.
Impact of different preprocessing pipelines on several metrics.
"""

import os
import scanpy as sc
import mito as mt
import matplotlib
import matplotlib.pyplot as plt
import plotting_utils as plu
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'data', 'tune')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')
path_results = os.path.join(path_main, 'results', 'others', 'Supp')


##


# Load tune results
df, metrics, options = mt.ut.format_tuning(path_data)
df.loc[df['pp_method']=='mito_preprocessing', 'pp_method'] = 'MiTo'

# Groups
groupings = ['pp_method', 'af_confident_detection', 'min_n_confidently_detected', 'bin_method', 'metric']
metric_annot = {
    'Mutation Quality' : ['n_dbSNP', 'n_REDIdb', 'transitions_vs_transversions_ratio'],
    'Association with GBC' : ['freq_lineage_biased_muts', 'AUPRC', 'ARI', 'NMI'],                               
    'Tree structure' : ['corr'],
    'Connectedness' : ['density', 'transitivity', 'average_path_length', 'average_degree', 'proportion_largest_component'],
    'Variation' : ['genomes_redundancy', 'median_n_vars_per_cell'],                                                           
    'Yield' : ['n_GBC_groups', 'n_cells', 'n_vars']                                                                
}  
weights = {
    'Mutation Quality': .1,
    'Association with GBC': .4,
    'Tree structure' : .1,
    'Connectedness' : .0,
    'Variation' : 0,
    'Yield' : .4
}
df_ranked = mt.ut.rank_items(df, groupings, metrics, weights, metric_annot)


##


# Colors
pp_method_colors = plu.create_palette(df, 'pp_method', sc.pl.palettes.vega_10)


# Paras
plu.set_rcParams()

# Figs
fig, ax = plt.subplots(figsize=(4,4))
plu.scatter(
    df_ranked, x='Association with GBC score', y='Yield score', by='pp_method', 
    categorical_cmap=pp_method_colors, size=10, ax=ax
)
plu.format_ax(ax=ax, xlabel='Association with GBC score', ylabel='Yield score')
plu.add_legend(
    ax=ax, colors=pp_method_colors, 
    loc='upper left', 
    bbox_to_anchor=(0,1), 
    label='Preprocessing\nmethod'
)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'Supp_fig_preprocessing_overall_scatter.pdf'))


##


fig, axs = plt.subplots(3,3,figsize=(8,6),sharex=True)

metrics_to_plot = [
    'n_cells', 'n_vars', 'median_n_vars_per_cell',
    'n_dbSNP', 'n_REDIdb', 'transitions_vs_transversions_ratio',
    'ARI', 'NMI', 'AUPRC'
]
metric_names = {'transitions_vs_transversions_ratio':'Transitions /\nTransversion', 
                 'median_n_vars_per_cell':'Median n_vars per cell'}

order = df_ranked.groupby('pp_method')['ARI'].median().sort_values().index
for i,metric in enumerate(metrics_to_plot):
    ax = axs.ravel()[i]
    df_feat = df.groupby(groupings, dropna=False)[metric].median().reset_index()
    name = metric_names[metric] if metric in metric_names else metric
    plu.box(df_feat, 'pp_method', metric, color='white', ax=ax, x_order=order)
    plu.format_ax(ax=ax, title=name, xlabel='', rotx=90, reduced_spines=True)

fig.subplots_adjust(bottom=.2, hspace=.4)
fig.savefig(os.path.join(path_figures, 'Supp_fig_preprocessing_metrics.pdf'))


##

