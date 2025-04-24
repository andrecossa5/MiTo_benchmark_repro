"""
Plasticity of gene expression programs (i.e., pySCENINC regulons) in PT-lung longitudinal dataset.
"""

import os
import pandas as pd
import numpy as np
import scanpy as sc
import pickle
import mito as mt
import matplotlib
import matplotlib.pyplot as plt
import plotting_utils as plu
from sklearn.preprocessing import scale
from mito.pp.filters import moran_I
from plotting_utils.plotting_base import _reorder
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_expression = os.path.join(path_main, 'data', 'longitudinal')
path_mt = os.path.join(path_main, 'results', 'others', 'Fig4', 'longitudinal')
path_others = os.path.join(path_main, 'results', 'others', 'Fig4', 'fitness')
path_figures = os.path.join(path_main, 'results', 'figures', 'Fig4')

# Read data
afm = sc.read(os.path.join(path_mt, 'afm_filtered.h5ad'))
expr = sc.read(os.path.join(path_expression, 'expression.h5ad'))
regulons = pd.read_csv(os.path.join(path_expression, 'regulons.csv'), index_col=0)
regulons['gene_set'] = regulons['gene_set'].map(lambda x: x.split(', '))
with open(os.path.join(path_mt, 'annotated_tree.pickle'), 'rb') as f:
    tree = pickle.load(f)

# Score regulons in single cells
scores = np.zeros((expr.shape[0], regulons.shape[0]))
for i,tf in enumerate(regulons.index):
    genes_set = regulons.loc[tf,'gene_set']
    if np.sum([ gene in expr.var_names for gene in genes_set ])>= len(genes_set) * .9:
        sc.tl.score_genes(expr, gene_list=genes_set)
        scores[:,i] = expr.obs['score'].values

# Filter out too weakly expressed regulons, and rescale them
test = np.sum(scores, axis=0)!=0
scores = pd.DataFrame(
    scale(scores[:,test]), 
    index=expr.obs_names, 
    columns=regulons.index[test]
)

# Test plasticity: Moran's I statistics 
W = 1/(1+afm.obsp['distances'].toarray())
I = []
P = [] 
for tf in scores.columns:
    x = scores.loc[afm.obs_names, tf].values
    i, p = moran_I(W, x)
    I.append(i)
    P.append(p)

stats = pd.DataFrame({'stat':I, 'pval':P, 'regulon':scores.columns})
stats = stats.sort_values('stat', ascending=False)
stats.describe()


##


# Visualization

# 1. Regulons plasticity (Moran's I) and cell state heamap
plu.set_rcParams({'figure.dpi':150})

fig, axs = plt.subplots(2,1,figsize=(6, 2.5), sharex=True, height_ratios=[1,2])

df_plot = scores.join(expr.obs[['cell_state']]).groupby('cell_state').median()
df_plot = df_plot.loc[['TGFb-EMT', 'IFN', 'Hypoxia', 'Proliferation', 'OXPHOS']]
order_cols = plu.order_from_index(df_plot)

stats['regulon'] = pd.Categorical(stats['regulon'], categories=order_cols)
plu.bar(stats, x='regulon', y='stat', ax=axs[0], color='white')
plu.format_ax(ax=axs[0], ylabel='Moran\'s I', reduced_spines=True)
plu.plot_heatmap(df_plot[order_cols], palette='viridis', 
                 label='Score', xlabel='pySCENIC regulons', cb=False,
                 ax=axs[1], y_names_size=10, x_names=True, x_names_size=3.5)
plu.add_cbar(df_plot.values.flatten(), ax=axs[1], palette='viridis', 
             label='Score', layout=( (1.025,.0,.01,1), 'right', 'vertical' ))

fig.subplots_adjust(top=.9, bottom=.3, left=.18, right=.85)
fig.savefig(os.path.join(path_figures, 'regulons_heatmap.pdf'))


##


# 2. Example of plastic and memory regulons
tree.cell_meta = tree.cell_meta.join(scores)
stats.head()


##


fig, axs = plt.subplots(2,4,figsize=(8,4.7))

TFs = [ 
    *stats.head(4)['regulon'].to_list(), 
    *stats.tail(4)['regulon'].to_list()
]
for i,tf in enumerate(TFs):
    ax = axs.ravel()[i]
    mt.pl.plot_tree(tree, ax=ax, 
                    continuous_cmaps={tf:'viridis'},
                    features=[tf],
                    orient=90,
                    vmin=-1.2, vmax=1.3,
                    colorstrip_width=25
                    )
    ax.set(title=tf)

plu.add_cbar(
    scores.values.flatten(), 
    ax=ax, palette='viridis', 
    vmin=-1.2, vmax=1.3,
    label='Score', 
    layout=( (1-.27,.05,.22,.02), 'bottom', 'horizontal' )
)
    
fig.tight_layout()
plu.save_best_pdf_quality(
    fig, 
    figsize=(8,4.7), 
    path=path_figures, 
    name='plasticity_regulons_trees.pdf'
)


##