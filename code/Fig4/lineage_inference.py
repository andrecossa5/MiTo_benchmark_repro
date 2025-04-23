"""
Prepare PT-lung couple for joint analysis. Visualize final tree.
"""

import os
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import mito as mt
import plotting_utils as plu
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_figures = os.path.join(path_main, 'results', 'figures', 'Fig4')


##


# 1. Prep afm with only cell-state annotated cells ------------------------------------#

# path_data = os.path.join(path_main, 'data', 'longitudinal')

# Read MT-SNVs and expression data
# afm = sc.read(os.path.join(path_data, 'afm_unfiltered.h5ad'))
# adata = sc.read(os.path.join(path_data, 'expression.h5ad'))
# 
# # Subset for common cells
# cells = list(set(adata.obs_names) & set(afm.obs_names))
# afm = afm[cells,:].copy()
# adata = adata[cells,:].copy()
# 
# # Reorder meta
# adata.obs['GBC'] = afm.obs['GBC']
# afm.obs['cell_state'] = adata.obs['cell_state']
# adata.obs = adata.obs[['GBC', 'sample', 'nUMIs', 'mito_perc', 
#            'detected_genes', 'doublet_score', 'leiden_0.11', 'cell_state']]
# afm.obs = afm.obs[['GBC', 'sample', 'nUMIs', 'mito_perc', 
#            'detected_genes', 'doublet_score', 
#            'mean_site_coverage', 'median_target_site_coverage', 
#            'median_untarget_site_coverage', 'frac_target_site_covered',
#            'cell_state']]
# 
# # Write out
# afm.write(os.path.join(path_data, 'afm_unfiltered.h5ad'))
# adata.write(os.path.join(path_data, 'expression.h5ad'))


##


# 2. Visualize lineage inference results ------------------------------------#

# Paths
path_data = os.path.join(path_main, 'results', 'others', 'Fig4', 'longitudinal')
path_figures = os.path.join(path_main, 'results', 'figures', 'Fig4')

# Read data
afm = sc.read(os.path.join(path_data, 'afm_filtered.h5ad'))
tree_metrics = pd.read_csv(os.path.join(path_data, 'tree_metrics.csv'), index_col=0)
with open(os.path.join(path_data, 'annotated_tree.pickle'), 'rb') as f:
    tree = pickle.load(f)


##


# Viz

# Params
plu.set_rcParams()
plt.rcParams.update({'figure.dpi': 150.0})  # Display with macosx

# fig, ax = plt.subplots(figsize=(3.5,3.5))
# mt.pl.heatmap_distances(afm, tree=tree, ax=ax)
# fig.tight_layout()
# plt.show()
# 
# fig, ax = plt.subplots(figsize=(3.5,3.5))
# mt.pl.heatmap_variants(afm, tree=tree, ax=ax, kwargs={'x_names_size':3}, cmap='afmhot_r')
# fig.tight_layout()
# plt.show()


##


fig, ax = plt.subplots(figsize=(6.5,5.5))
cmaps = {
    'MiTo clone':plu.create_palette(tree.cell_meta, 'MiTo clone', sc.pl.palettes.default_102),
    'cell_state':plu.create_palette(tree.cell_meta, 'cell_state', col_list=plu.darjeeling),
    'sample' : {'MDA_PT':'#00928E', 'MDA_lung':'#DC4C0D'}
}
tree.cell_meta['MiTo clone'] = pd.Categorical(tree.cell_meta['MiTo clone'])
mt.pl.plot_tree(
    tree, ax=ax, 
    features=['MiTo clone', 'cell_state', 'sample'], 
    orient='down',
    categorical_cmaps=cmaps, 
    colorstrip_width=6,
    label_offset=5,
    label_size=10,
    feature_internal_nodes='support',
    cmap_internal_nodes='Spectral_r',
    internal_node_subset=tree.cell_meta['lca'].value_counts().index,
    internal_node_kwargs={'markersize':5},
    vmin_internal_nodes=.5, vmax_internal_nodes=.8
)
tree_stats = mt.tl.get_internal_node_stats(tree)
plu.add_cbar(tree_stats['support'], palette='Spectral_r', ax=ax, vmin=.5, vmax=.8,
             layout=( (.75,1,.2,.015), 'top', 'horizontal' ), label='Support', ticks_size=5)

fig.subplots_adjust(left=.2, right=.9, bottom=.1, top=.85)
fig.savefig(os.path.join(path_figures, 'tree.png'), dpi=1000)
 

##