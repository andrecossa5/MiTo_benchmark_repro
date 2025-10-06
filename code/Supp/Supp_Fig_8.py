"""
Supp Fig 8
Plot maegatk/MiTo with and without MiTo filter.
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


# 2. Plot -------------------------- 

x_order = ['MDA_clones', 'MDA_lung', 'MDA_PT']
by_order = ['MiTo-MQuad', 'MiTo-MiTo', 'maegatk-MQuad', 'maegatk-MiTo']
colors = ['#EE8383', '#d62728', '#D99BEC', '#aa40fc']
cmap = dict(zip(by_order, colors))


##

fig, axs = plt.subplots(1,2,figsize=(7,3.5))

plu.bar(
    df.loc[df['combo'].isin(by_order)].set_index('job_id'), 
    x='sample', y='ARI', 
    by='combo',
    x_order=x_order,
    by_order=by_order,
    categorical_cmap=cmap,
    ax=axs[0]
)
plu.format_ax(ax=axs[0], xlabel='', ylabel='ARI', reduced_spines=True, rotx=90)
plu.bar(
    df.loc[df['combo'].isin(by_order)].set_index('job_id'), 
    x='sample', y='NMI', 
    by='combo',
    x_order=x_order,
    by_order=by_order,
    categorical_cmap=cmap,
    ax=axs[1]
)
plu.format_ax(ax=axs[1], xlabel='', ylabel='NMI', reduced_spines=True, rotx=90)
plu.add_legend(cmap, label='Preprocessing-\nMT-SNVs filtering\ncombination', ax=axs[1])

fig.subplots_adjust(top=.90, bottom=.3, right=.7, left=.1, wspace=.35)
fig.savefig(os.path.join(path_figures, 'Supp_Fig_8.pdf'))


##
