"""
Supp Fig 2
Bench cell_filter, preliminary parameter tuning results.
"""

import os
import pandas as pd
import mito as mt
import matplotlib
import matplotlib.pyplot as plt
import plotting_utils as plu
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'data', 'bench', 'tune_filters_preliminary')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')


# Set visualization params
plu.set_rcParams({'figure.dpi':350})


##


# 1. Extract data from TUNE output with preliminary parameter tuning ------------------------------------#

# Format data
L = []
for folder,_,files in os.walk(path_data):
    if any([ x.startswith('all') for x in files]):
        df,metrics,options = mt.ut.format_tuning(folder)
        L.append(df)
df = pd.concat(L)

# Define options and metrics
varying_options = (df[options].nunique()).loc[lambda x:x>1].index.to_list()
varying_options
metrics_of_interest = ['ARI', 'NMI', 'corr', 'AUPRC', 'n_cells', 'n_vars', 'n_GBC_groups']

##

# Get the ordered list of all unique combinations of varying options
# with filter 2 (i.e., best cell filter)
df_combos = (
    df
    .query('cell_filter=="filter2"')
    .assign(combo = lambda x: x[varying_options].astype(str).agg('-'.join, axis=1))
    [['combo', 'sample']+metrics_of_interest]
)
x_order = df_combos.groupby('combo')['ARI'].median().sort_values(ascending=False).index
samples = ['MDA_clones', 'MDA_PT', 'MDA_lung']


##


# 2. Supp Fig 2. key metrics for each grid search combo ------------------------------------#

df = df.query('cell_filter=="filter2"')
sample_colors = plu.create_palette(df, 'sample', 'tab10')
n_unique = df_combos['combo'].unique().size
n_trials = df_combos.shape[0]
x_label = f'Unique combinations of (n=5) cell and MT-SNV filtering options\n(n unique combinations={n_unique}, n total trials={n_trials})'
metrics_of_interest = ['ARI', 'AUPRC', 'n_cells', 'n_GBC_groups']

fig, axs = plt.subplots(4,1,figsize=(10,5))

for i,metric in enumerate(metrics_of_interest):

    ax = axs.ravel()[i]
    plu.strip(
        df_combos, 
        x='combo', y=metric, 
        categorical_cmap=sample_colors,
        x_order=x_order, 
        by_order=samples, 
        by='sample', 
        ax=ax, 
        size=3
    )
    ax.axhline(df.query('sample=="MDA_clones"')[metric].median(), linewidth=.8, c=sample_colors["MDA_clones"])
    ax.axhline(df.query('sample=="MDA_PT"')[metric].median(), linewidth=.8, c=sample_colors["MDA_PT"])
    ax.axhline(df.query('sample=="MDA_lung"')[metric].median(), linewidth=.8, c=sample_colors["MDA_lung"])

    xticks = [''for _ in range(len(x_order))]
    plu.format_ax(
        ax=ax, xticks=xticks, ylabel=metric, 
        xlabel='' if metric != 'n_GBC_groups' else x_label,
        reduced_spines=True
    )
    ax.grid(axis='x', linewidth=.2)

    if metric == 'ARI':
        plu.add_legend(
            colors=sample_colors, 
            ax=ax, 
            label='sample', 
            loc='lower left', bbox_to_anchor=(0,0),
            ncols=3
        )

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'Supp_Fig_2a.pdf'))


##


# 3. Supp Fig 2. cell filters ------------------------------------#

fig = plt.figure(figsize=(10,4))

df = pd.concat(L)
cmap = plu.create_palette(df, 'cell_filter', palette='Blues')
metrics_of_interest = ['ARI', 'NMI', 'corr', 'AUPRC', 'n_cells', 'n_vars', 'n_GBC_groups']

for i,y in enumerate(metrics_of_interest):
    ax = fig.add_subplot(2,4,i+1)
    plu.bar(df, x='sample', y=y, by='cell_filter', categorical_cmap=cmap, ax=ax)
    plu.format_ax(
        reduced_spines=True, 
        xlabel='', 
        ylabel=y,
        xticks=['','',''] if i<4 else None,
        rotx=90,
        ax=ax
    )

plu.add_legend(cmap, label='Cell filter', ax=ax)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'Supp_Fig_2b.pdf'))


##