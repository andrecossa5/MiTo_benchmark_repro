"""
Bench cell_filter, preliminary. (Supp. Fig. )
"""

import os
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
path_data = os.path.join(path_main, 'data', 'bench', 'tune_filters_preliminary')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')


##


# Format data
L = []
for folder,_,files in os.walk(path_data):
    if any([ x.startswith('all') for x in files]):
        df,metrics,options = mt.ut.format_tuning(folder)
        L.append(df)
df = pd.concat(L)


# Define options and metrics
varying_options = (df[options].nunique()).loc[lambda x:x>1].index.to_list()
metrics_of_interest = ['ARI', 'NMI', 'corr', 'AUPRC', 'n_cells', 'n_vars', 'n_GBC_groups']


##


# Plot
plu.set_rcParams()
plt.rcParams.update({'figure.dpi':350})


##


# Combos
df_combos = (
    df.query('cell_filter=="filter2"')
    .assign(combo = lambda x: x[varying_options].astype(str).agg('-'.join, axis=1))
    [['combo', 'sample']+metrics_of_interest]
)
x_order = df_combos.groupby('combo')['ARI'].median().sort_values(ascending=False).index
samples = ['MDA_clones', 'MDA_PT', 'MDA_lung']

##


fig, axs = plt.subplots(4,1,figsize=(10,5))

df = df.query('cell_filter=="filter2"')
sample_colors = plu.create_palette(df, 'sample', 'tab10')

plu.strip(df_combos, x='combo', y='ARI', categorical_cmap=sample_colors,
          x_order=x_order, by_order=samples, by='sample', ax=axs[0], size=3)
axs[0].axhline(df.query('sample=="MDA_clones"')['ARI'].median(), linewidth=.8, c=sample_colors["MDA_clones"])
axs[0].axhline(df.query('sample=="MDA_PT"')['ARI'].median(), linewidth=.8, c=sample_colors["MDA_PT"])
axs[0].axhline(df.query('sample=="MDA_lung"')['ARI'].median(), linewidth=.8, c=sample_colors["MDA_lung"])
axs[0].set(ylim=(0,1))
plu.format_ax(ax=axs[0], xticks=[''for _ in range(len(x_order))], 
              xlabel='', ylabel='ARI', reduced_spines=True)
axs[0].grid(axis='x', linewidth=.2)
plu.add_legend(
    colors=sample_colors, 
    ax=axs[0], 
    label='sample', 
    loc='lower left', bbox_to_anchor=(0,0),
    ncols=3
)

plu.strip(df_combos, x='combo', y='AUPRC', categorical_cmap=sample_colors,
          x_order=x_order, by_order=samples, by='sample', ax=axs[1], size=3)
axs[1].axhline(df.query('sample=="MDA_clones"')['AUPRC'].median(), linewidth=.8, c=sample_colors["MDA_clones"])
axs[1].axhline(df.query('sample=="MDA_PT"')['AUPRC'].median(), linewidth=.8, c=sample_colors["MDA_PT"])
axs[1].axhline(df.query('sample=="MDA_lung"')['AUPRC'].median(), linewidth=.8, c=sample_colors["MDA_lung"])
axs[1].set(ylim=(0,1))
plu.format_ax(ax=axs[1], xticks=[''for _ in range(len(x_order))], 
              xlabel='', ylabel='AUPRC', reduced_spines=True)
axs[1].grid(axis='x', linewidth=.2)

plu.strip(df_combos, x='combo', y='n_cells', categorical_cmap=sample_colors,
          x_order=x_order, by_order=samples, by='sample', ax=axs[2], size=3)
axs[2].axhline(df.query('sample=="MDA_clones"')['n_cells'].median(), linewidth=.8, c=sample_colors["MDA_clones"])
axs[2].axhline(df.query('sample=="MDA_PT"')['n_cells'].median(), linewidth=.8, c=sample_colors["MDA_PT"])
axs[2].axhline(df.query('sample=="MDA_lung"')['n_cells'].median(), linewidth=.8, c=sample_colors["MDA_lung"])
plu.format_ax(ax=axs[2], xticks=[''for _ in range(len(x_order))], 
              xlabel='', ylabel='n_cells', reduced_spines=True)
axs[2].grid(axis='x', linewidth=.2)

plu.strip(df_combos, x='combo', y='n_GBC_groups', 
          x_order=x_order, by_order=samples, by='sample', ax=axs[3], size=3)
axs[3].axhline(df.query('sample=="MDA_clones"')['n_GBC_groups'].median(), linewidth=.8, c=sample_colors["MDA_clones"])
axs[3].axhline(df.query('sample=="MDA_PT"')['n_GBC_groups'].median(), linewidth=.8, c=sample_colors["MDA_PT"])
axs[3].axhline(df.query('sample=="MDA_lung"')['n_GBC_groups'].median(), linewidth=.8, c=sample_colors["MDA_lung"])
# axs[1].set(ylim=(200,3000))
n_unique = df_combos['combo'].unique().size
n_trials = df_combos.shape[0]
plu.format_ax(ax=axs[3], xticks=[''for _ in range(len(x_order))], 
              xlabel=f'Unique combinations of (n=5) cell and MT-SNV filtering options\n(n unique combinations={n_unique}, n total trials={n_trials})',
              ylabel='n_GBC_groups', reduced_spines=True)
axs[3].grid(axis='x', linewidth=.2)

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'Supp_Fig_2_up.pdf'))


##


# Cell filters, global performance tested over other 4 major variant filter varying options.
fig = plt.figure(figsize=(10,4))

df = pd.concat(L)
cmap = plu.create_palette(df, 'cell_filter', palette='Blues')
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
fig.savefig(os.path.join(path_figures, 'Supp_Fig_2_down.pdf'))


##