"""
Supp Fig ... binarization strategies comparison.
Impact of different binarization strategies on several metrics.
"""

import os
import numpy as np
import pandas as pd
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


# 1. Jobs combos. All jobs included (only filter2). -------------------------- 

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
bin_method_colors = plu.create_palette(df, 'bin_method', 'Reds')


# Paras
plu.set_rcParams()

# Figs
fig, ax = plt.subplots(figsize=(4,4))
plu.scatter(
    df_ranked, x='Association with GBC score', y='Yield score', by='bin_method', 
    categorical_cmap=bin_method_colors, size=10, ax=ax
)
plu.format_ax(ax=ax, xlabel='Association with GBC score', ylabel='Yield score')
plu.add_legend(
    ax=ax, colors=bin_method_colors, 
    loc='upper left', 
    bbox_to_anchor=(0,1), 
    label='Genotyping \nmethod'
)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'Supp_fig_binarization_overall_scatter.pdf'))


##


fig, axs = plt.subplots(3,3,figsize=(8,6),sharex=True)

metrics_to_plot = [
    'n_cells', 'n_vars', 'median_n_vars_per_cell', 'corr', 
    'mean_CI', 'ARI', 'NMI', 'AUPRC', 'freq_lineage_biased_muts'

]
metric_names = {'transitions_vs_transversions_ratio':'Transitions /\nTransversion', 
                 'median_n_vars_per_cell':'Median n_vars\nper cell',
                 'corr':'Tree-char dists\ncorrelation',
                 'freq_lineage_biased_muts' : '% biased MT-SNVs'
                }

x_order = df_ranked.groupby('pp_method')['ARI'].median().sort_values().index
by_order = ['vanilla', 'MiTo_smooth', 'MiTo']

for i,metric in enumerate(metrics_to_plot):
    ax = axs.ravel()[i]
    df_feat = df.groupby(groupings, dropna=False)[metric].median().reset_index()
    name = metric_names[metric] if metric in metric_names else metric
    plu.box(
        df_feat, 'pp_method', metric, by='bin_method', ax=ax,
        categorical_cmap=bin_method_colors, x_order=x_order, by_order=by_order
    )
    plu.format_ax(ax=ax, title=name, xlabel='', rotx=90, reduced_spines=True)

fig.subplots_adjust(bottom=.2, hspace=.4)
fig.savefig(os.path.join(path_figures, 'Supp_fig_bin_metrics_metrics.pdf'))


##


# 2. Single jobs. Only maegatk, filter2, MiTo. -------------------------- 
path_data = os.path.join(path_main, 'results', 'others', 'Fig2')

# Read selected jobs
L = []
for x in os.listdir(path_data):
    if x.endswith('filtered_jobs.csv'):
        L.append(pd.read_csv(os.path.join(path_data, x), index_col=0))
df = pd.concat(L)

# Filter only: pp_method == maegatk, cell_filter==filter2, filtering==MiTo
df = (
    df.loc[
        (df['cell_filter']=='filter2') & \
        (df['pp_method'].isin(['maegatk'])) & \
        (df['filtering']=='MiTo') 
    ]
)

##


fig, axs = plt.subplots(1,3,figsize=(6,3.5),sharey=True)

metric = 'NMI'
bin_method_colors = {
    'vanilla': (0.9882352941176471, 0.732072279892349, 0.6299269511726259),
    'MiTo': (0.7925720876585928, 0.09328719723183392, 0.11298731257208766)
}

for i,sample in enumerate(['MDA_clones', 'MDA_PT', 'MDA_lung']):
    ax = axs[i]
    plu.bar(df.query('sample==@sample'), 'af_confident_detection', metric, 
            categorical_cmap=bin_method_colors, by='bin_method', by_order=['vanilla', 'MiTo'], ax=ax, alpha=1)
    plu.format_ax(ax=ax, title=sample, reduced_spines=True, xlabel='Min confident AF', ylabel=metric)

fig.tight_layout()
fig.savefig(os.path.join(path_figures, f'Supp_fig_binarization_{metric}_selected_jobs.pdf'))


##
