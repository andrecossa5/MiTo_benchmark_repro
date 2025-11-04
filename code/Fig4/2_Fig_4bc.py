"""
Fig 4:
- Harmonize PT-lung trimodal dataset. 
- Visualization of:
  1. Joint PT-lung gene expression space (Fig 4b)
  2. Joint PT-lung cell phylogeny (Fig 4c)
"""

import os
import anndata
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
path_afms = os.path.join(path_main, 'data', 'general', 'AFMs')
path_longitudinal = os.path.join(path_main, 'data', 'longitudinal')
path_lineage_inference = os.path.join(path_main, 'results', 'others', 'Fig4', 'lineage_inference')
path_figures = os.path.join(path_main, 'results', 'figures', 'Fig4')

# Set visualization params
plu.set_rcParams({'figure.dpi':350})


##


# 1. Prep merged, unfiltered AFM matrix for the PT-lung dataset -------------------------- #

# Merge AFMs
afm_PT = sc.read(os.path.join(path_afms, 'MDA_PT', 'afm_unfiltered.h5ad'))
afm_lung = sc.read(os.path.join(path_afms, 'MDA_lung', 'afm_unfiltered.h5ad'))
afm = anndata.concat([afm_PT, afm_lung])

# Fix .var and .uns
afm.var['pos'] = afm.var_names.map(lambda x: x.split('_')[0]).astype('int')
afm.var['ref'] = afm.var_names.map(lambda x: x.split('_')[1].split('>')[0])
afm.var['alt'] = afm.var_names.map(lambda x: x.split('_')[1].split('>')[1])
afm.uns['pp_method'] = 'maegatk'
afm.uns['scLT_system'] = 'MAESTER'
afm


##


# 2. Subset GEX and MT data for the same set of cells ------------------------------------#

# Read GEX data
adata = sc.read(os.path.join(path_longitudinal, 'expression.h5ad'))

# Subset for common cells
cells = list(set(adata.obs_names) & set(afm.obs_names))
afm = afm[cells,:].copy()
adata = adata[cells,:].copy()

# Re-format meta

# Add missing columns, reciprocally
adata.obs['GBC'] = afm.obs['GBC']
afm.obs['cell_state'] = adata.obs['cell_state']

# Reorder columns, keep only relevant ones
adata.obs = adata.obs[
    ['GBC', 'sample', 'nUMIs', 
     'mito_perc', 'leiden_0.11', 'cell_state']
]
afm.obs = afm.obs[
    ['GBC', 'sample', 'nUMIs', 'mito_perc', 
     'mean_site_coverage', 'median_target_site_coverage', 
     'median_untarget_site_coverage', 'frac_target_site_covered',
     'cell_state']
]

# Write out harmonized AnnData objects
afm.write(os.path.join(path_longitudinal, 'afm_unfiltered.h5ad'))
adata.write(os.path.join(path_longitudinal, 'expression.h5ad'))


##

# 3. Lineage inference ------------------------------------#

"""
nf-MiTo workflow goes here.
See nf-MiTo INFER workflow at https://github.com/andrecossa5/nf-MiTo.
Results are provided in the <source data>/results/others/Fig4/lineage_inference folder.
"""

##


# 4. Fig 4b. Gene expression space visualization ------------------------------------#

# Read filtered AFM and expression data
afm = sc.read(os.path.join(path_lineage_inference, 'afm_filtered.h5ad'))
adata = sc.read(os.path.join(path_longitudinal, 'expression.h5ad'))

# Subset GEX data ONLY for cells present in the final filtered AFM 
adata = adata[afm.obs_names,:].copy()

# Set colors
cmaps = {
    'cell_state':plu.create_palette(adata.obs, 'cell_state', col_list=plu.darjeeling),
    'sample' : {'MDA_PT':'#00928E', 'MDA_lung':'#DC4C0D'}
}

# Plot embeddings
fig, axs = plt.subplots(1,2,figsize=(4.3,2.2))
for i,feat in enumerate(['sample', 'cell_state']):
    sc.pl.umap(adata, color=feat, ax=axs[i], show=False, 
               palette=cmaps[feat], legend_loc='on data', frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'Fig_2b.pdf'))


##


# 5. Fig 4c. Annotated single-cell phylogeny ------------------------------------#

# Load phylogeny
tree_metrics = pd.read_csv(os.path.join(path_lineage_inference, 'tree_metrics.csv'), index_col=0)
with open(os.path.join(path_lineage_inference, 'annotated_tree.pickle'), 'rb') as f:
    tree = pickle.load(f)

# Plot tree
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
tree_stats = mt.ut.get_internal_node_stats(tree)
plu.add_cbar(tree_stats['support'], palette='Spectral_r', ax=ax, vmin=.5, vmax=.8,
             layout=( (.75,1,.2,.015), 'top', 'horizontal' ), label='Support', ticks_size=5)

fig.subplots_adjust(left=.2, right=.9, bottom=.1, top=.85)
plu.save_best_pdf_quality(
    fig, 
    figsize=(6.5,5.5),
    path=path_figures,
    name='Fig_4c.pdf'
)


##