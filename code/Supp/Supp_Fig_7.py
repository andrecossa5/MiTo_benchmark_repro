"""
Supp Fig 7
--pp-method AND --filtering method combinations benchmark.
"""

import os
import pandas as pd
import mito as mt
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import plotting_utils as plu
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'data', 'bench', 'tune_pp_method')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')


# Set visualization params
plu.set_rcParams({'figure.dpi':350})


##


# 1. Extract data from TUNE output with --pp-method and --filtering tuning ------------------------------------#

L = []
for folder,_,files in os.walk(path_data):
    if any([ x.startswith('all') for x in files]):
        df,metrics,options = mt.ut.format_tuning(folder)
        L.append(df)
df = pd.concat(L)
df.loc[df['pp_method'] == 'mito_preprocessing', 'pp_method'] = 'MiTo'
df['combo'] = df[['pp_method', 'filtering']].astype(str).agg('-'.join, axis=1)

# Define options and metrics
varying_options = (df[options].nunique()).loc[lambda x:x>1].index.to_list()
metrics_of_interest = ['ARI', 'NMI', 'corr', 'AUPRC', 'n_cells', 'n_vars', 'n_GBC_groups']
metrics_of_interest += [
    'freq_lineage_biased_muts', 'median_n_vars_per_cell', 
    'transitions_vs_transversions_ratio', 'n_dbSNP', 'n_REDIdb'
]


##


# 2. Plot ------------------------------------#

# Pre-processing method - filtering combo
fig =  plt.figure(figsize=(14,6.5))

x_order = ['MDA_clones', 'MDA_lung', 'MDA_PT']
by_order = ['samtools-baseline', 'freebayes-baseline', 'cellsnp-lite-MQuad', 'MiTo-MiTo', 'maegatk-MiTo']
cmap = plu.create_palette(df, 'combo', order=by_order, col_list=sc.pl.palettes.vega_10_scanpy)
test = (df['combo'].isin(['maegatk-MQuad', 'MiTo-MQuad']))

for i,metric in enumerate(metrics_of_interest):
    ax = fig.add_subplot(2,6,i+1)
    plu.bar(
        df.loc[~test].set_index('job_id'), 
        x='sample', y=metric, 
        by='combo',
        x_order=x_order,
        by_order=by_order,
        categorical_cmap=cmap,
        ax=ax
    )
    plu.format_ax(ax=ax, 
                  xticks=x_order if i>=6 else [ '' for _ in range(3) ], 
                  rotx=90,
                  xlabel='', ylabel=metric, reduced_spines=True)
    if i==5:
        plu.add_legend(cmap, label='Preprocessing-\nMT-SNVs filtering\ncombination', 
                       ax=ax, bbox_to_anchor=(1,0.3), loc='upper left')

fig.subplots_adjust(top=.90, bottom=.25, right=.82, left=.05, wspace=.8)
fig.savefig(os.path.join(path_figures, 'Supp_fig_7.pdf'))


##