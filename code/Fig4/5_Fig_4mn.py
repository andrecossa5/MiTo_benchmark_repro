"""
Fig 4:
- Infer the plasticity of gene expression programs (i.e., pySCENINC regulons)
- Visualization of:
  1. pySCENIC regulons, their Moran's I, and cell state bias (Fig 4m)
  2. Examples of memory and plasticity regulons (Fig 4n)
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
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_GEX = os.path.join(path_main, 'data', 'longitudinal')
path_lineage_inference = os.path.join(path_main, 'results', 'others', 'Fig4', 'lineage_inference')
path_figures = os.path.join(path_main, 'results', 'figures', 'Fig4')

# Set visualization params
plu.set_rcParams({'figure.dpi':350})


##


# 1. Score pySCENIC regulons, compute Morans'I statistic ------------------------------------#

# Read data
afm = sc.read(os.path.join(path_lineage_inference, 'afm_filtered.h5ad'))
adata = sc.read(os.path.join(path_GEX, 'expression.h5ad'))
with open(os.path.join(path_lineage_inference, 'annotated_tree.pickle'), 'rb') as f:
    tree = pickle.load(f)
regulons = pd.read_csv(os.path.join(path_GEX, 'regulons.csv'), index_col=0)
regulons['gene_set'] = regulons['gene_set'].map(lambda x: x.split(', '))

# Score regulons in single cells (i.e., sc.tl.score_genes)
scores = np.zeros((adata.shape[0], regulons.shape[0]))
for i,tf in enumerate(regulons.index):
    genes_set = regulons.loc[tf,'gene_set']
    genes_present = [ gene in adata.var_names for gene in genes_set ]
    test = np.sum(genes_present)>= len(genes_set) * .9
    if test:
        sc.tl.score_genes(adata, gene_list=genes_set)
        scores[:,i] = adata.obs['score'].values

# Filter out too weakly expressed regulons, and z-score their expression score
test = np.sum(scores, axis=0)!=0
scores = pd.DataFrame(scale(scores[:,test]), index=adata.obs_names, columns=regulons.index[test])

# Compute Moran's I statistics 
W = 1/(1+afm.obsp['distances'].toarray())
I = []
P = [] 
for tf in scores.columns:
    x = scores.loc[afm.obs_names, tf].values
    i, p = mt.pp.filters.moran_I(W, x)
    I.append(i)
    P.append(p)

# Store Moran's I stat values
stats = pd.DataFrame({'stat':I, 'pval':P, 'regulon':scores.columns})
stats = stats.sort_values('stat', ascending=False)
stats.describe()


##


# 2. Fig 4m. pySCENIC regulons and their Morans'I stats ---------------------#

# Plot
fig, axs = plt.subplots(2,1,figsize=(6, 2.5), sharex=True, height_ratios=[1,2])

df_plot = scores.join(adata.obs[['cell_state']]).groupby('cell_state').median()
df_plot = df_plot.loc[['TGFb-EMT', 'IFN', 'Hypoxia', 'Proliferation', 'OXPHOS']]
order_cols = plu.order_from_index(df_plot)

stats['regulon'] = pd.Categorical(stats['regulon'], categories=order_cols)
plu.bar(stats, x='regulon', y='stat', ax=axs[0], color='white')
plu.format_ax(ax=axs[0], ylabel='Moran\'s I', reduced_spines=True)
plu.plot_heatmap(
    df_plot[order_cols], palette='viridis', 
    label='Score', xlabel='pySCENIC regulons', cb=False,
    ax=axs[1], y_names_size=10, x_names=True, x_names_size=3.5
)
plu.add_cbar(df_plot.values.flatten(), ax=axs[1], palette='viridis', 
             label='Score', layout=( (1.025,.0,.01,1), 'right', 'vertical' ))

fig.subplots_adjust(top=.9, bottom=.3, left=.18, right=.85)
fig.savefig(os.path.join(path_figures, 'Fig_4m.pdf'))


##


# 3. Fig 4n. Examples of memory and plasticity regulons ---------------------#

# Join pySCENIC scores to cell phylogeny cell metadata
tree.cell_meta = tree.cell_meta.join(scores)

# Plot
fig, axs = plt.subplots(2,4,figsize=(8,4.7))

TFs = [ 
    *stats.head(4)['regulon'].to_list(), 
    *stats.tail(4)['regulon'].to_list()
]
for i,tf in enumerate(TFs):
    ax = axs.ravel()[i]
    mt.pl.plot_tree(
        tree, 
        ax=ax, 
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
    name='Fig_4n.pdf'
)


##