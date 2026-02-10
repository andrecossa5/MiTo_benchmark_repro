"""
Supp Fig 19:
- Supplementary Gene Expression data for longitudinal dataset.
"""

import os
import numpy as np
import scanpy as sc
import mito as mt
import plotting_utils as plu
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_longitudinal = os.path.join(path_main, 'data', 'longitudinal')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')


# Set visualization params
plu.set_rcParams({'figure.dpi':350})


##


# Read filtered AFM and expression data
adata = sc.read(os.path.join(path_longitudinal, 'expression.h5ad'))


##


# Supp Fig 19a. Dotplot cell_state markers ------------------------------------#

# Marker genes
group = 'cell_state'
min_n_cells = adata.obs[group].value_counts().min() * .5
test = np.sum(adata.layers['raw'].toarray()>0, axis=0) >= min_n_cells
adata = adata[:,test].copy()
sc.tl.rank_genes_groups(adata, groupby=group, method='wilcoxon', pts=True)

# Filter final marker genes
res_de = mt.ut.format_rank_genes_groups(adata)
order_groups = mt.ut.order_groups(adata, groupby=group, obsm_key='X_pca', n_dims=20)
res_de_filtered = mt.ut.format_rank_genes_groups(adata, filter_genes=True)
top_markers = mt.ut.get_top_markers(res_de_filtered, order_groups=order_groups, ntop=3)
df_plot = res_de.query('gene in @top_markers')

# Dotplot
fig, ax = plt.subplots(figsize=(6,3))
plu.dotplot(df_plot, 'gene', 'group', 
            order_x=top_markers, order_y=order_groups,
            color='log2FC', size='pct_group', ax=ax, vmin=-5, vmax=5)
plu.format_ax(ax=ax, rotx=90, xlabel='', ylabel='Clusters')
ax.get_legend().set_bbox_to_anchor((1,1.5))
ax.margins(x=0.1, y=0.2)
fig.subplots_adjust(top=.65, bottom=.35, left=.25, right=.7)
fig.savefig(os.path.join(path_figures, 'Fig_19a.pdf'))


##


# Supp Fig 19b. Cell state abundance ------------------------------------#

fig, ax = plt.subplots(figsize=(3.5,1.3))
cmap = plu.create_palette(adata.obs, 'cell_state', col_list=plu.darjeeling)
plu.bb_plot(adata.obs, 'sample', 'cell_state', categorical_cmap=cmap, ax=ax)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'Fig_19b.pdf'))


##


# Supp Fig 19c. UMAPs technical covariates ------------------------------------#

# Add detected genes to metadata
adata.obs['detected_genes'] = (adata.X>0).sum(1).A1

# Plot embeddings
fig, axs = plt.subplots(1,3,figsize=(6,2))
for i,feat in enumerate(['nUMIs', 'mito_perc', 'detected_genes']):
    sc.pl.umap(adata, color=feat, ax=axs[i], show=False, frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'Fig_19c.pdf'))


##