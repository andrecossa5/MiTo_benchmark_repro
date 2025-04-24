"""
Supp Fig ... Full MiTo and nf-MiTo visualization. 
Shocase plot gallery.
"""

import os
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import mito as mt
import matplotlib
import matplotlib.pyplot as plt
import plotting_utils as plu
from plotting_utils.plotting_base import _reorder
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'data', 'bench', 'clonal_inference')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')
path_results = os.path.join(path_main, 'results', 'others', 'Supp')

# Select sample and job
sample = 'MDA_clones'
job_id = '92bbd42de6'
path_data = os.path.join(path_data, sample, job_id)

# Load: afm, tree, metrics, phylocorr
afm = sc.read(os.path.join(path_data, 'afm.h5ad'))
with open(os.path.join(path_data, 'annotated_tree.pickle'), 'rb') as f:
    tree = pickle.load(f)
metrics = pd.read_csv(os.path.join(path_data, 'tree_metrics.csv'), index_col=0).set_index('metric')
phylocorr = pd.read_csv(os.path.join(path_data, 'phylocorr.csv'), index_col=0)

# Load colors
path_colors = os.path.join(path_main, 'data', 'general')
with open(os.path.join(path_colors, 'clones_colors_sc.pickle'), 'rb') as f:
    colors = pickle.load(f)


##


# Viz
plu.set_rcParams()

# UMAP
mt.pp.reduce_dimensions(afm)

# Main UMAP
fig, ax = plt.subplots(figsize=(4,4))
mt.pl.draw_embedding(afm, feature='GBC', ax=ax, size=250, categorical_cmap=colors)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, f'Supp_full_viz_{sample}_{job_id}_embeddings.pdf'))


##


# Select clones and variants
mapping = {
    'GTTCGAAGACTCCAGGGA' : '2732_G>A',
    'CGGAGCAGCCGCAGCGGG' : '2659_C>T',
    'CACTGAAGCTGTTCGGCG' : '11072_T>C'
}

fig, axs = plt.subplots(2,3,figsize=(9,6.5))

for i,clone in enumerate(mapping):

    mut = mapping[clone]
    afm.obs['cat'] = np.where(afm.obs['GBC']==clone, clone, 'unassigned')
    _cmap = { clone : colors[clone], 'unassigned' : 'lightgrey' }
    mt.pl.draw_embedding(afm, feature='cat', categorical_cmap=_cmap, ax=axs[0,i], size=150)
    mt.pl.draw_embedding(afm, feature=mut, ax=axs[1,i], continuous_cmap='mako_r', 
                         size=150, legend=False, kwargs={'colorbar_loc':None})
    axs[0,i].set(title=clone[:5])
    axs[1,i].set(title=mut)

fig.subplots_adjust(bottom=.2)
fig.savefig(os.path.join(path_figures, f'Supp_full_viz_{sample}_{job_id}_embeddings.pdf'))


##


# Trees
cmaps = { 
    'GBC' : colors,
    'MiTo clone' : plu.create_palette(tree.cell_meta, 'MiTo clone', 
                                      col_list=sc.pl.palettes.vega_10_scanpy, add_na=True)
}
tree.cell_meta['MiTo clone'][tree.cell_meta['MiTo clone'].isna()] = 'unassigned'
clonal_nodes = tree.cell_meta['lca'].unique()[1:]
tree_stats = mt.tl.get_internal_node_stats(tree)

# Fig
fig, axs = plt.subplots(1,2,figsize=(8,4.2))

mt.pl.plot_tree(
    tree, 
    categorical_cmaps=cmaps, 
    ax=axs[0], 
    colorstrip_width=5,
    internal_node_subset=clonal_nodes,
    feature_internal_nodes='support',
    internal_node_kwargs={'markersize':7.5}
)
plu.add_cbar(
    tree_stats['support'], palette='Spectral_r', ax=axs[0], 
    vmin=.5, vmax=.8, label='Support', layout='outside'
)
corr = metrics.loc['corr_distances', 'value']
supp_mut_clades = metrics.loc['median_support_mut_clades', 'value']
title = f'Clonal nodes support: {supp_mut_clades:.2f}\nTree-character distance corr: {corr:.2f}'
plu.format_ax(ax=axs[0], title=title)

##

mt.pl.plot_tree(
    tree, 
    features=['GBC', 'MiTo clone'], 
    categorical_cmaps=cmaps, 
    ax=axs[1], 
    colorstrip_width=5,
    internal_node_subset=clonal_nodes,
    show_internal=True,
    internal_node_kwargs={'markersize':7.5, 'c':'darkorange'}
)
ari = metrics.loc['ARI', 'value']
nmi = metrics.loc['NMI', 'value']
title = f'ARI: {ari:.2f}, NMI: {nmi:.2f}'
plu.format_ax(ax=axs[1], title=title)

fig.tight_layout()
fig.savefig(os.path.join(path_figures, f'Supp_full_viz_{sample}_{job_id}_trees.pdf'))


##

fig, ax = plt.subplots(figsize=(4.5,4))

df_corr = phylocorr.pivot(index='Var1', columns='Var2', values='Z')
order = _reorder(df_corr.values)
ax.imshow(df_corr.iloc[order, order], cmap='magma')
plu.format_ax(title='PATH phylo-corr', ax=ax, xticks=[], yticks=[], xlabel='GBC', ylabel='GBC')
plu.add_cbar(x=df_corr.values.flatten(), label='phylocorr', ax=ax, palette='magma')
fig.tight_layout()
fig.savefig(os.path.join(path_figures, f'Supp_full_viz_{sample}_{job_id}_phylocorr.pdf'))


##


fig, axs = plt.subplots(1,2,figsize=(11,5))

model = mt.tl.MiToTreeAnnotator(tree)
model.get_T()
model.get_M()
model.extract_mut_order()

mt.pl.plot_tree(
    tree, 
    features=['GBC', 'MiTo clone'], 
    categorical_cmaps=cmaps, 
    characters=model.ordered_muts,
    layer='raw',
    orient='down',
    ax=axs[0], 
    colorstrip_width=5,
    internal_node_subset=clonal_nodes,
    feature_internal_nodes='support',
    labels=True,
    internal_node_kwargs={'markersize':7.5}
)
plu.add_cbar(
    afm.X.A.flatten(), palette='mako', ax=axs[0], 
    vmin=.01, vmax=.1, label='Allelic Frequency', layout='outside'
)

##

mt.pl.plot_tree(
    tree,
    features=['GBC', 'MiTo clone'],  
    categorical_cmaps=cmaps, 
    characters=model.ordered_muts,
    layer='transformed',
    orient='down',
    ax=axs[1], 
    colorstrip_width=5,
    internal_node_subset=clonal_nodes,
    feature_internal_nodes='support',
    labels=False,
    internal_node_kwargs={'markersize':7.5}
)
plu.add_legend(
    colors={'MUT':'r','WT':'b' },
    ax=axs[1], label='Genotype', loc='center left', bbox_to_anchor=(1,.5)
)

fig.subplots_adjust(wspace=.4, left=.15, right=.85)
fig.savefig(os.path.join(path_figures, f'Supp_full_viz_{sample}_{job_id}_mut_trees.png'), dpi=1000)


##