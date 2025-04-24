"""
Supp Fig ... . Accuracy-vs-molecular detection metrics.
Impact of filtering options controlling sensitivity of molecular detection and 
clonal inference accuracy.
"""

import os
import numpy as np
import pandas as pd
import mito as mt
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import plotting_utils as plu
from itertools import product
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

# Define columns, and names to plot
columns_ = [
    'af_confident_detection', 'min_mean_AD_in_positives', 
    'min_AD', 'median_n_vars_per_cell'
]
col_names_ = {
    'af_confident_detection' : 'AF confident detection',
    'min_mean_AD_in_positives' : 'Mean n ALT UMIs in +cells', 
    'min_AD' : 'n ALT UMIs\nfor genotyping',
    'median_n_vars_per_cell' : 'Median n MT-SNVs per cell'
}
GBC_cols = ['NMI', 'ARI', 'AUPRC']


##


# Viz
plu.set_rcParams()

# Fig
fig, axs = plt.subplots(3,4,figsize=(9,7.5), sharey=True, gridspec_kw={'width_ratios':[5,3,2,4]})

df_plot = (
    df[['pp_method', 'bin_method']+columns_+GBC_cols]
    .dropna()
    .groupby(['pp_method', 'bin_method']+columns_)
    [GBC_cols].mean()
    .reset_index()
)

for i,(metric,col) in enumerate(list(product(GBC_cols, columns_))):
    ax = axs.ravel()[i]
    if col == 'median_n_vars_per_cell':
        sns.regplot(data=df_plot, x=col, y=metric, ax=ax, scatter=False)
        ax.plot(df_plot[col], df_plot[metric], 
                color='lightgrey', marker='o', linestyle='', 
                markersize=5, markeredgecolor='k')
    else:
        plu.box(df_plot, col, metric, ax=ax, 
                color='white', x_order=sorted(df_plot[col].unique()))
        
    name = col_names_[col] if col in col_names_ else col
    name = name if i in [8,9,10,11] else ''
    metric = metric if i in [0,4,8] else ''
    plu.format_ax(ax=ax, xlabel=name, ylabel=metric, reduced_spines=True)

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'Supp_fig_accuracy_vs_molecular_detection.pdf'))


##