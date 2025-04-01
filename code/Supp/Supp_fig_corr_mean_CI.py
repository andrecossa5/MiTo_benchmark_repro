"""
Supp Fig ... criteria to inform model selection.
Correlation among tree structural metrics and GBC-accuracy prediction metrics.
"""

import os
import numpy as np
import mito as mt
import matplotlib
import matplotlib.pyplot as plt
import plotting_utils as plu
import seaborn as sns
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


# Params
plu.set_rcParams()


##


## Viz
fig, axs = plt.subplots(1,4,figsize=(13,3.5))

metrics = ['freq_lineage_biased_muts', 'AUPRC', 'ARI', 'NMI']

for i,metric in enumerate(metrics):
    ax = axs.ravel()[i]
    sns.regplot(data=df_ranked, x='mean_CI', y=metric, ax=ax, scatter=False)
    ax.plot(df_ranked['mean_CI'], df_ranked[metric], color='lightgrey', 
            marker='o', linestyle='', markersize=5, markeredgecolor='k')
    corr = np.corrcoef(df_ranked['mean_CI'].values, df_ranked[metric].values)[0,1]
    ax.set(title=f'Pearson\'s r: {corr:.2f}')

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'mean_CI_and_other_GBC.pdf'))


##


fig, axs = plt.subplots(1,4,figsize=(13,3.5))

metrics = ['freq_lineage_biased_muts', 'AUPRC', 'ARI', 'NMI']

for i,metric in enumerate(metrics):
    ax = axs.ravel()[i]
    sns.regplot(data=df_ranked, x='corr', y=metric, ax=ax, scatter=False)
    ax.plot(df_ranked['corr'], df_ranked[metric], color='lightgrey', 
            marker='o', linestyle='', markersize=5, markeredgecolor='k')
    corr = np.corrcoef(df_ranked['corr'].values, df_ranked[metric].values)[0,1]
    ax.set(title=f'Pearson\'s r: {corr:.2f}')

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'corr_and_other_GBC.pdf'))


##