"""
Basic gene expression analysis for PT-lung longitudinal dataset.
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
from anndata import AnnData
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'data', 'longitudinal')
path_figures = os.path.join(path_main, 'results', 'figures', 'Fig4')

# Read expression data
adata = sc.read(os.path.join(path_data, 'QC.h5ad'))
adata.layers['raw'] = adata.X.copy()


##


# Exploratory: spot bad cell populations
n_cells = adata.shape[0]
robust_genes = np.sum(adata.X.A>0, axis=0)> (.01 * n_cells)
test = adata.var_names.str.startswith('MT-') | \
       adata.var_names.str.startswith('RPL') | \
       adata.var_names.str.startswith('RPS')

adata = adata[:,(robust_genes) & (~test)].copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
sc.tl.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, metric='euclidean')
sc.tl.umap(adata)

res = np.linspace(0.1,2,10)
for r in res:
    sc.tl.leiden(adata, resolution=r, key_added=f'leiden_{r:.2f}')

plu.set_rcParams({'figure.dpi':100})
sc.pl.umap(adata, color=['sample'])
sc.pl.umap(adata, color=adata.obs.columns[adata.obs.columns.str.startswith('leiden')])

# Remove bad clusters: leiden_0.10, 1 , 2
adata.obs.groupby('leiden_0.10')[['nUMIs', 'detected_genes']].describe().T
adata = adata[~adata.obs['leiden_0.10'].isin(["1","2"])].copy()

# Restore anndata
adata = AnnData(X=adata.layers['raw'], obs=adata.obs.iloc[:,:5], var=pd.DataFrame(index=adata.var_names))
adata.layers['raw'] = adata.X.copy()


##


# Single-cell workflow

# Gene filtering
n_cells = adata.shape[0]
robust_genes = np.sum(adata.X.A>0, axis=0)> (.01 * n_cells)
test = adata.var_names.str.startswith('MT-') | \
       adata.var_names.str.startswith('RPL') | \
       adata.var_names.str.startswith('RPS')

adata = adata[:,(robust_genes) & (~test)].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")

# Remove cc genes from highly_variable_genes
ki67_index = np.where(adata.var_names=='MKI67')[0]
corr = np.corrcoef(adata.X.A.T)

# sns.kdeplot(corr[:,ki67_index])
# plt.show()

adata.var['cell_cycle'] = corr[:,ki67_index].flatten()>.5
adata.var['highly_variable'] = (adata.var['highly_variable']) & (~adata.var['cell_cycle'])

sc.tl.pca(adata, n_comps=50)
sc.pl.pca_variance_ratio(adata)
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=20, metric='euclidean')
sc.tl.diffmap(adata, n_comps=10)
sc.pp.neighbors(adata, n_neighbors=50, use_rep="X_diffmap",metric='euclidean')
sc.tl.umap(adata)
sc.tl.draw_graph(adata)

# Viz embeddings
sc.pl.umap(adata, color=['sample'])
sc.pl.draw_graph(adata, color=['sample'])
sc.pl.diffmap(adata, color=['sample'])

# Leiden clustering
res = np.linspace(0.05,.15,10)
for r in res:
    sc.tl.leiden(adata, resolution=r, key_added=f'leiden_{r:.2f}')

plu.set_rcParams({'figure.dpi':100})
sc.pl.umap(adata, color=adata.obs.columns[adata.obs.columns.str.startswith('leiden')])
sc.pl.draw_graph(adata, color=adata.obs.columns[adata.obs.columns.str.startswith('leiden')])

# Choose best solution: silhouette score
# from sklearn.metrics import silhouette_score
# silhouettes = []
# for col in adata.obs.columns[adata.obs.columns.str.startswith('leiden')]:
#     silhouettes.append(silhouette_score(adata.obsp['distances'].A, adata.obs[col].values))
# plt.plot(adata.obs.columns[adata.obs.columns.str.startswith('leiden')], silhouettes, 'ko')
# plt.show()

# Choose leiden solution
chosen = 'leiden_0.11'

# Viz chosen solution
sc.pl.umap(adata, color=chosen)
sc.pl.draw_graph(adata, color=chosen)

# Polish adata
adata.obs = adata.obs.drop(
    columns=adata.obs.columns[ 
        (adata.obs.columns.str.startswith('leiden')) & (adata.obs.columns != chosen) 
    ]
)
uns = {}
for x in adata.uns:
    if not x.endswith('colors'):
        uns[x] = adata.uns[x]
adata.uns = uns


##


# Cell state annotation

# Viz embeddings
sc.pl.umap(adata, color=chosen)
sc.pl.umap(adata, color='nUMIs')
sc.pl.umap(adata, color='detected_genes')
sc.pl.umap(adata, color='mito_perc')
sc.pl.diffmap(adata, color='sample')
sc.pl.umap(adata, color='sample')

##

# Annotate cell states from unsupervised clustering

# QC metrics
adata.obs.groupby(chosen).median()

# Markers
adata.uns['log1p']['base'] = np.exp(1)
min_n_cells = adata.obs[chosen].value_counts().min() * .5
test = (~adata.var_names.str.contains('.', regex=False)) & \
       (~adata.var_names.str.contains('LIN')) & \
       (np.sum(adata.layers['raw'].A>0, axis=0) >= min_n_cells)
adata = adata[:,test].copy()
sc.tl.rank_genes_groups(adata, groupby=chosen, method='wilcoxon', pts=True)

# Format and filter markers genes
res_de = mt.ut.format_rank_genes_groups(adata)
order_groups = mt.ut.order_groups(adata, groupby=chosen, obsm_key='X_pca', n_dims=20)
res_de_filtered = mt.ut.format_rank_genes_groups(adata, filter_genes=True)
top_markers = mt.ut.get_top_markers(res_de_filtered, order_groups=order_groups, ntop=3)
df_plot = res_de.query('gene in @top_markers')

# Dotplot top markers
fig, ax = plt.subplots(figsize=(6,3.5))
plu.dotplot(df_plot, 'gene', 'group', 
            order_x=top_markers, order_y=order_groups,
            color='log2FC', size='pct_group', ax=ax, vmin=-5, vmax=5)
plu.format_ax(ax=ax, rotx=90, xlabel='', ylabel='Clusters')
ax.get_legend().set_bbox_to_anchor((1,1.2))
fig.tight_layout()
plt.show()

# Spot individual markers and associated biological processes
cluster = "4" # Example, leiden_0.11: cluster 4
res_de_filtered.query('group==@cluster')
res_gsea = mt.ut.run_GSEA(
       res_de.query('group==@cluster').set_index('gene')['log2FC'], 
    max_size_set=1000, 
    collection='MSigDB_Hallmark_2020'
)
res_gsea['Term'] = res_gsea['Term'].map(lambda x: x.replace('MSigDB_Hallmark_2020__', ''))
res_gsea

##

# Cell state annotation
mapping_cell_states = {
    '0' : 'OXPHOS',
    '1' : 'TGFb-EMT',
    '2' : 'Proliferation',
    '3' : 'TGFb-EMT',
    '4' : 'IFN',
    '5' : 'Hypoxia',
    '6' : np.nan,    # Poor quality, undefined biological processes
    '7' : np.nan,    
    '8' : np.nan    
}
adata.obs['cell_state'] = adata.obs[chosen].map(mapping_cell_states)

# Remove undefined/poor-quality cells
adata = adata[~adata.obs['cell_state'].isna(),:].copy()

# Remove strange group of cells from TGFb-EMT cluster
test = (adata.obsm['X_umap'][:,0] < 5) & \
       (adata.obsm['X_umap'][:,1] > 5) & \
       (adata.obs['cell_state'] == 'TGFb-EMT')
adata.obs.loc[test,'cell_state'] = np.nan
adata = adata[~adata.obs['cell_state'].isna(),:].copy()

# Re-calculate markers
chosen = 'cell_state'
adata.uns['log1p']['base'] = np.exp(1)
min_n_cells = adata.obs[chosen].value_counts().min() * .5
test = (~adata.var_names.str.contains('.', regex=False)) & \
       (~adata.var_names.str.contains('LIN')) & \
       (np.sum(adata.layers['raw'].A>0, axis=0) >= min_n_cells)
adata = adata[:,test].copy()
sc.tl.rank_genes_groups(adata, groupby=chosen, method='wilcoxon', pts=True)

# Format and filter markers genes
res_de = mt.ut.format_rank_genes_groups(adata)
order_groups = mt.ut.order_groups(adata, groupby=chosen, obsm_key='X_pca', n_dims=20)
res_de_filtered = mt.ut.format_rank_genes_groups(adata, filter_genes=True)
top_markers = mt.ut.get_top_markers(res_de_filtered, order_groups=order_groups, ntop=3)
df_plot = res_de.query('gene in @top_markers')

# Write out
adata.write(os.path.join(os.path.join(path_data, 'expression.h5ad')))


##


# Visualization
adata = sc.read(os.path.join(os.path.join(path_data, 'expression.h5ad')))

# Params
plu.set_rcParams({'figure.dpi':150})

# Embeddings
fig, axs = plt.subplots(1,2,figsize=(5,2.5))
for i,feat in enumerate(['sample', 'cell_state']):
       sc.pl.umap(adata, color=feat, ax=axs[i], show=False, legend_loc='on data', frameon=False)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'expression_umap.pdf'))

##

fig, axs = plt.subplots(1,3,figsize=(8,2.5))
for i,feat in enumerate(['nUMIs', 'mito_perc', 'detected_genes']):
       sc.pl.umap(adata, color=feat, ax=axs[i], show=False, frameon=False)
fig.tight_layout()
plt.show()

##

# Abundance cell states
fig, ax = plt.subplots(figsize=(3.5,1.3))
plu.bb_plot(adata.obs, 'sample', 'cell_state', ax=ax)
fig.tight_layout()
plt.show()

# Marker genes
fig, ax = plt.subplots(figsize=(6,3))
plu.dotplot(df_plot, 'gene', 'group', 
            order_x=top_markers, order_y=order_groups,
            color='log2FC', size='pct_group', ax=ax, vmin=-5, vmax=5)
plu.format_ax(ax=ax, rotx=90, xlabel='', ylabel='Clusters')
ax.get_legend().set_bbox_to_anchor((1,1.5))
ax.margins(x=0.1, y=0.2)
fig.subplots_adjust(top=.65, bottom=.35, left=.25, right=.7)
plt.show()


##







