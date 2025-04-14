"""
Clonal bias and evolutionary coupling.
"""

import os
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import plotting_utils as plu
import matplotlib
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from sklearn.metrics import pairwise_distances
from moscot.problems.time import LineageProblem
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


# 2. Visualize clonal bias as cell state occupancy ------------------------------------#

# Paths
path_expression = os.path.join(path_main, 'data', 'longitudinal')
path_mt = os.path.join(path_main, 'results', 'others', 'Fig4', 'longitudinal')
path_figures = os.path.join(path_main, 'results', 'figures', 'Fig4')

# Read data
afm = sc.read(os.path.join(path_mt, 'afm_filtered.h5ad'))
tree_metrics = pd.read_csv(os.path.join(path_mt, 'tree_metrics.csv'), index_col=0)
with open(os.path.join(path_mt, 'annotated_tree.pickle'), 'rb') as f:
    tree = pickle.load(f)


##


# Viz

# Params
plu.set_rcParams()
plt.rcParams.update({'figure.dpi': 150.0})  # Display with macosx

# Frequencies
fig, axs = plt.subplots(2,1,figsize=(5,2.5))

cell_states = tree.cell_meta['cell_state'].value_counts().index.values.to_list()
plu.plot_heatmap(
    pd.crosstab(tree.cell_meta['cell_state'], tree.cell_meta['MiTo clone'], normalize=1).loc[cell_states], 
    palette='Blues', ax=axs[0], cluster_cols=True, cluster_rows=False, x_names=False, xlabel='MiTo clones',
    label='Frequency'
)
plu.plot_heatmap(
    pd.crosstab(tree.cell_meta['cell_state'], tree.cell_meta['GBC'], normalize=1).loc[cell_states], 
    palette='Blues', ax=axs[1], cluster_cols=True, cluster_rows=False, x_names=False, xlabel='GBC clones',
    label='Frequency'
)

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'clonal_freqs.pdf'))


## 


# 3. Moslin transition rates ------------------------------------#

# Paths
path_expression = os.path.join(path_main, 'data', 'longitudinal')
path_mt = os.path.join(path_main, 'results', 'others', 'Fig4', 'longitudinal')
path_others = os.path.join(path_main, 'results', 'others', 'Fig4')
path_figures = os.path.join(path_main, 'results', 'figures', 'Fig4')

# Read data
afm = sc.read(os.path.join(path_mt, 'afm_filtered.h5ad'))
tree_metrics = pd.read_csv(os.path.join(path_mt, 'tree_metrics.csv'), index_col=0)
with open(os.path.join(path_mt, 'annotated_tree.pickle'), 'rb') as f:
    tree = pickle.load(f)
adata = sc.read(os.path.join(path_expression, 'expression.h5ad'))

# Subset adata
adata = adata[afm.obs_names].copy()

##

# Run moslin: MT

# Add MT-distances
adata.obsp['mt_distances'] = afm.obsp['distances']
# Format time
adata.obs['time'] = adata.obs['sample'].map({'MDA_PT':1, 'MDA_lung':2})

# Moslin
lp = LineageProblem(adata)
lp = lp.prepare(
    time_key="time",
    joint_attr="X_pca",
    lineage_attr={"attr": "obsp", "key": "mt_distances", "cost": "custom"},
)
lp = lp.solve(alpha=0.99, epsilon=1e-3, tau_a=0.99)
cell_states = adata.obs['cell_state'].value_counts().index.values.to_list()
cell_transition = lp.cell_transition(
    source=1,
    target=2,
    source_groups={"cell_state": cell_states},
    target_groups={"cell_state": cell_states},
    forward=True,
    key_added="lp_transitions",
)
cell_transition.to_csv(os.path.join(path_others, 'moslin_cell_transitions_MT.csv'))

##

# Run moslin: GBC

# Remove previous results
del adata.uns['moscot_results']

# Add GBC-distances
adata.obs['GBC'] = adata.obs['GBC'].astype(str)
adata.obs['value'] = 1
GBC_bin = (
    adata.obs.reset_index()
    .pivot_table(index='index', values='value', columns='GBC', fill_value=0)
)
D = pairwise_distances(GBC_bin.values, metric='jaccard')
adata.obsp['gbc_distances'] = csr_matrix(D)

# Moslin
lp = LineageProblem(adata)
lp = lp.prepare(
    time_key="time",
    joint_attr="X_pca",
    lineage_attr={"attr": "obsp", "key": "gbc_distances", "cost": "custom"},
)
lp = lp.solve(alpha=0.99, epsilon=1e-3, tau_a=0.99)
cell_states = adata.obs['cell_state'].value_counts().index.values.to_list()
cell_transition = lp.cell_transition(
    source=1,
    target=2,
    source_groups={"cell_state": cell_states},
    target_groups={"cell_state": cell_states},
    forward=True,
    key_added="lp_transitions",
)
cell_transition.to_csv(os.path.join(path_others, 'moslin_cell_transitions_gbc.csv'))


##


# Viz estimated cell transition
mt_tr = pd.read_csv(os.path.join(path_others, 'moslin_cell_transitions_mt.csv'), index_col=0)
gbc_tr = pd.read_csv(os.path.join(path_others, 'moslin_cell_transitions_gbc.csv'), index_col=0)

# Params
plu.set_rcParams()
plt.rcParams.update({'figure.dpi': 150.0})  # Display with macosx

# Frequencies
fig, axs = plt.subplots(1,2,figsize=(4.5,2.5))

plu.plot_heatmap(mt_tr, palette='inferno', ax=axs[0], x_names=True, y_names=True, label='p')
axs[0].set(title='MT-based')
plu.plot_heatmap(gbc_tr, palette='inferno', ax=axs[1], x_names=True, y_names=False, label='p')
axs[1].set(title='GBC-based')

fig.tight_layout()

fig.savefig(os.path.join(path_figures, 'cell_transitions.pdf'))


##