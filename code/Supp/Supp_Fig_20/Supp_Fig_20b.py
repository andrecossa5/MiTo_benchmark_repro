"""
Supp Fig 20
- Re-analysis of Weng et. al, 2024 dataset: mouse batch 1.
- Visualize metrics from TUNE of RedeeM preprocessing and MT-SNVs filtering.
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


# Set visualization params
plu.set_rcParams({'figure.dpi':120})


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'data', 'redeem_mouse')
path_results = os.path.join(path_main, 'results', 'others', 'Supp20')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')


##


# Format tuning
df, _, _ = mt.ut.format_tuning(os.path.join(path_results, 'TUNE'))

# Select best MT-SNVs spaces (highest ARI, with at least 75 cells per sample)
# L = []
# grouped = (
#     df.query('filtering=="MiTo" and n_cells>=75')
#     .groupby('sample')
#     [['filtering', 'trim', 'job_id', 'ARI', 'n_vars', 'n_cells']]
# )
# for name, group in grouped:
#     L.append(group.sort_values('ARI', ascending=False).head(1).assign(sample=name))
# df_selected = pd.concat(L).reset_index(drop=True)[['job_id', 'sample', 'filtering', 'trim']]
# df_selected.to_csv(os.path.join(path_results, 'selected_jobs_Supp20cd.csv'), index=False)


##


# Visualize results
fig, axs = plt.subplots(2,3,figsize=(7,4.5), sharex=True)

colors = plu.create_palette(df, 'trim', palette=plu.darjeeling)
x_order = ['baseline', 'RedeemR', 'MiTo']

plu.box(df, 'filtering', 'n_cells', ax=axs[0,0], by='trim', categorical_cmap=colors, x_order=x_order)
plu.format_ax(ax=axs[0,0], xlabel='', ylabel='n cells', reduced_spines=True, rotx=90)

plu.box(df, 'filtering', 'n_vars', ax=axs[0,1], by='trim', categorical_cmap=colors, x_order=x_order)
plu.format_ax(ax=axs[0,1], xlabel='', ylabel='n MT-SNVs (total)', reduced_spines=True, rotx=90)

plu.box(df, 'filtering', 'mean_n_vars_per_cell', ax=axs[0,2], by='trim', categorical_cmap=colors, x_order=x_order)
plu.format_ax(ax=axs[0,2], xlabel='', ylabel='n MT-SNVs per cell', reduced_spines=True, rotx=90)
plu.add_legend(colors, label='Genotypes', ax=axs[0,2])

plu.box(df, 'filtering', 'ARI', ax=axs[1,0], by='trim', categorical_cmap=colors, x_order=x_order)
plu.format_ax(ax=axs[1,0], xlabel='', reduced_spines=True, rotx=90)

plu.box(df, 'filtering', 'freq_lineage_biased_muts', ax=axs[1,1], by='trim', categorical_cmap=colors, x_order=x_order)
plu.format_ax(ax=axs[1,1], xlabel='', ylabel='% lineage-biased MT-SNVs', reduced_spines=True, rotx=90)

plu.box(df, 'filtering', 'corr', ax=axs[1,2], by='trim', categorical_cmap=colors, x_order=x_order)
plu.format_ax(ax=axs[1,2], xlabel='', ylabel='Distance correlation', reduced_spines=True, rotx=90)

fig.subplots_adjust(right=.8, left=.2, top=.9, bottom=.2, wspace=.8)
fig.savefig(os.path.join(path_figures, 'Supp_20b.pdf'), dpi=300)


##