"""
Pyscenic regulons PT-lung longitudinal dataset.
"""

import os
import pandas as pd
import numpy as np
import scanpy as sc
import mito as mt
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import plotting_utils as plu
from sklearn.preprocessing import scale
from plotting_utils.plotting_base import _reorder
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'data', 'longitudinal')
path_figures = os.path.join(path_main, 'results', 'figures', 'Fig4')


##


# Read expression data
expr = sc.read(os.path.join(path_data, 'expression.h5ad'))

# Get pyscenic regulons 
regulons = pd.read_csv(os.path.join(path_data, 'regulons.csv'), index_col=0)
regulons['gene_set'] = regulons['gene_set'].map(lambda x: x.split(', '))

# Score regulons
scores = np.zeros((expr.shape[0], regulons.shape[0]))
for i,tf in enumerate(regulons.index):
    genes_set = regulons.loc[tf,'gene_set']
    if np.sum([ gene in expr.var_names for gene in genes_set ])>= len(genes_set) * .9:
        sc.tl.score_genes(expr, gene_list=genes_set)
        scores[:,i] = expr.obs['score'].values

test = np.sum(scores, axis=0)!=0
df_plot = pd.DataFrame(
    scale(scores[:,test]), index=expr.obs_names, columns=regulons.index[test]
)
df_plot = df_plot.join(expr.obs[['cell_state']]).groupby('cell_state').median()


##


# Viz
plu.set_rcParams()

fig, ax = plt.subplots(figsize=(6, 1.5))

cell_states = expr.obs['cell_state'].value_counts().index
df_plot = df_plot.loc[cell_states]
order_cols = plu.order_from_index(df_plot)

plu.plot_heatmap(df_plot[order_cols], palette='afmhot_r', 
                 cluster_cols=True, cluster_rows=True, 
                 label='Score', xlabel='pySCENIC regulons',  
                 ax=ax, y_names_size=10, x_names=True, x_names_size=3.5)

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'regulons_heatmap.pdf'))


##

