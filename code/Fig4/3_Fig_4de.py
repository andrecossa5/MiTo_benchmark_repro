"""
Fig 4:
- Compute PT-lung cell state transition rates with Moslin
- Visualization of:
  1. MT-vs-GBC clonal biases (Fig 4d)
  2. MT-vs-GBC transition rates (Fig 4e)
"""

import os
import pickle
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
path_GEX = os.path.join(path_main, 'data', 'longitudinal')
path_lineage_inference = os.path.join(path_main, 'results', 'others', 'Fig4', 'lineage_inference')
path_transition_rates = os.path.join(path_main, 'results', 'others', 'Fig4', 'transition_rates')
path_figures = os.path.join(path_main, 'results', 'figures', 'Fig4')

# Set visualization params
plu.set_rcParams({'figure.dpi':350})


##


# 1. Fig 4d. Visualize clonal bias as cell state occupancy ------------------------------------#

# Read annotated cell phylogeny
with open(os.path.join(path_lineage_inference, 'annotated_tree.pickle'), 'rb') as f:
    tree = pickle.load(f)

# Plot
fig, axs = plt.subplots(2,1,figsize=(5,2.5))

# Order cell states
cell_states = tree.cell_meta['cell_state'].value_counts().index.values.to_list()

# MT
df_ = pd.crosstab(tree.cell_meta['cell_state'], tree.cell_meta['MiTo clone'], normalize=1)
df_ = df_.loc[cell_states]
plu.plot_heatmap(
    df_.loc[:,plu.order_from_index(df_)], 
    palette='Blues', ax=axs[0], x_names=False, xlabel='MiTo clones',
    label='Frequency'
)

# GBC
df_ = pd.crosstab(tree.cell_meta['cell_state'], tree.cell_meta['GBC'], normalize=1)
df_ = df_.loc[cell_states]
plu.plot_heatmap(
    df_.loc[:,plu.order_from_index(df_)], 
    palette='Blues', ax=axs[1], x_names=False, xlabel='GBC clones',
    label='Frequency'
)

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'Fig_4d.pdf'))


## 


# 2. Infer PT-lung cell state transition rates with Moslin ------------------------------------#

# Read filtered AFM and GEX data
afm = sc.read(os.path.join(path_lineage_inference, 'afm_filtered.h5ad'))
adata = sc.read(os.path.join(path_GEX, 'expression.h5ad'))

# Subset adata for ONLY cells present in the filtered AFM
adata = adata[afm.obs_names].copy()


##


# 2.1. Run moslin with input distances in MT space

# Add MT-distances to adata, format time-series column as numeric
adata.obsp['mt_distances'] = afm.obsp['distances']
adata.obs['time'] = adata.obs['sample'].map({'MDA_PT':1, 'MDA_lung':2})

# Moslin, MT
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
cell_transition.to_csv(os.path.join(path_transition_rates, 'moslin_cell_transitions_MT.csv'))


##


# 2.2. Run moslin with input distances in GBC space

# Remove previous results
del adata.uns['moscot_results']

# Format GBC column as string
adata.obs['GBC'] = adata.obs['GBC'].astype(str)

# Calculate distances in GBC space 
"""
Here we will consider these distances as binary distances between each cell pair: 
Same clone: 1, different clone: 0.
"""

GBC_bin = (
    adata.obs
    .assign(value=1)
    .reset_index()
    .pivot_table(index='index', values='value', columns='GBC', fill_value=0)
)
D = pairwise_distances(GBC_bin.values, metric='jaccard')
adata.obsp['gbc_distances'] = csr_matrix(D)

# Moslin, GBC
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
cell_transition.to_csv(os.path.join(path_transition_rates, 'moslin_cell_transitions_gbc.csv'))


##


# 3. Fig 4.e. Visualize transition matrices ------------------------------------#

mt_tr = pd.read_csv(os.path.join(path_transition_rates, 'moslin_cell_transitions_mt.csv'), index_col=0)
gbc_tr = pd.read_csv(os.path.join(path_transition_rates, 'moslin_cell_transitions_gbc.csv'), index_col=0)

# Visualize transition rate matrices
fig, axs = plt.subplots(1,2,figsize=(5,2.5))

plu.plot_heatmap(mt_tr, palette='Spectral_r', ax=axs[0], x_names=True, y_names=True, label='p')
axs[0].set(title='MT-based')
plu.plot_heatmap(gbc_tr, palette='Spectral_r', ax=axs[1], x_names=True, y_names=False, label='p')
axs[1].set(title='GBC-based')

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'Fig_4e.pdf'))


##
