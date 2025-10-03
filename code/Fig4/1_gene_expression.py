"""
Fig 4:
- Gene expression analysis for the PT-lung longitudinal dataset.
"""

import os
import pandas as pd
import numpy as np
import scanpy as sc
import mito as mt
import matplotlib
import matplotlib.pyplot as plt
import plotting_utils as plu
from anndata import AnnData
# matplotlib.use('macOSX')        # On macOS only 


##


# Set paths 
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'data', 'longitudinal')
path_results = os.path.join(path_main, 'results', 'others', 'Fig4')

# Set visualization params
plu.set_rcParams({'figure.dpi':100})

##


# Read expression data (after first, sample specific lenient QC, 
# #and matrix concatenation)
adata = sc.read(os.path.join(path_data, 'QC.h5ad'))
adata.layers['raw'] = adata.X.copy()


##


# 1. Refine cell QC: spot additional low-quality cell populations --------------------------------------- #

# Remove genes too lowly expressed or are MT/ribosomal
n_cells = adata.shape[0]
robust_genes = np.sum(adata.X.toarray()>0, axis=0)> (.01 * n_cells)
test = adata.var_names.str.startswith('MT-') | \
       adata.var_names.str.startswith('RPL') | \
       adata.var_names.str.startswith('RPS')
adata = adata[:,(robust_genes) & (~test)].copy()

# Classic scanpy pipeline for single-cell GEX data preprocessing
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
sc.tl.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, metric='euclidean')
sc.tl.umap(adata)

# Clustering at different resolutions
res = np.linspace(0.1,2,10)
for r in res:
    sc.tl.leiden(adata, resolution=r, key_added=f'leiden_{r:.2f}')

# Visualize embeddings (samples and leiden clusters)
sc.pl.umap(adata, color=['sample'])
sc.pl.umap(adata, color=adata.obs.columns[adata.obs.columns.str.startswith('leiden')])

# Assess 2 suspicious clusters, even at low resolution (i.e., 0.10)
adata.obs.groupby('leiden_0.10')[['nUMIs', 'detected_genes']].describe().T

# Remove these cell populations (low number of UMIs and genes detected)
adata = adata[~adata.obs['leiden_0.10'].isin(["1","2"])].copy()

# Restore Anndata to repeat the single-cell workflow again
adata = AnnData(
    X=adata.layers['raw'], 
    obs=adata.obs.iloc[:,:5], 
    var=pd.DataFrame(index=adata.var_names),
    layers={'raw':adata.X.copy()}
)


##


# 2. Unsupervised clustering --------------------------------------- #

# Remove genes too lowly expressed or are MT/ribosomal
n_cells = adata.shape[0]
robust_genes = np.sum(adata.X.toarray()>0, axis=0)> (.01 * n_cells)
test = adata.var_names.str.startswith('MT-') | \
       adata.var_names.str.startswith('RPL') | \
       adata.var_names.str.startswith('RPS')
adata = adata[:,(robust_genes) & (~test)].copy()

# Log-normalization, HVG selection
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")

# Assess cell cycle correlation on GEX data
ki67_index = np.where(adata.var_names=='MKI67')[0]
corr = np.corrcoef(adata.X.toarray().T)

# Visualize gene-gene correlation distribution (with MKI67 gene)
# import seaborn as sns
# sns.kdeplot(corr[:,ki67_index])
# plt.show()
 
# Exlclude MKI67-correlated genes from HVGs
adata.var['cell_cycle'] = corr[:,ki67_index].flatten()>.5
adata.var['highly_variable'] = (adata.var['highly_variable']) & (~adata.var['cell_cycle'])

# Dimensionality reduction (with diffmap-denoised kNN graph)
sc.tl.pca(adata, n_comps=50)
sc.pl.pca_variance_ratio(adata) # n_pcs chosen=20
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=20, metric='euclidean')
sc.tl.diffmap(adata, n_comps=10)
sc.pp.neighbors(adata, n_neighbors=50, use_rep="X_diffmap",metric='euclidean')
sc.tl.umap(adata)

# Visualize new cell embeddings
sc.pl.umap(adata, color=['sample', 'nUMIs', 'detected_genes', 'mito_perc'])

# Leiden clustering
res = np.linspace(0.05,.15,10)
for r in res:
    sc.tl.leiden(adata, resolution=r, key_added=f'leiden_{r:.2f}')
# Visualize leiden clusters
sc.pl.umap(adata, color=adata.obs.columns[adata.obs.columns.str.startswith('leiden')])

# Choose best solution with silhouette score (long)
# from sklearn.metrics import silhouette_score
# silhouettes = []
# for col in adata.obs.columns[adata.obs.columns.str.startswith('leiden')]:
#     silhouettes.append(silhouette_score(adata.obsp['distances'].toarray(), adata.obs[col].values))
# plt.plot(adata.obs.columns[adata.obs.columns.str.startswith('leiden')], silhouettes, 'ko')
# plt.show()

# Chosen resolution: 'leiden_0.11'
chosen = 'leiden_0.11'
sc.pl.umap(adata, color=chosen)

# Polish adata .obs and .uns
adata.obs = (
    adata.obs.loc[:,~adata.obs.columns.str.startswith('leiden')]
    .join(adata.obs[[chosen]])
)
uns = {}
for x in adata.uns:
    if not x.endswith('colors'):
        uns[x] = adata.uns[x]
adata.uns = uns

# Assess QC metrics in chosen clustering solution
adata.obs.groupby(chosen).describe().T

# Find cluster markers
adata.uns['log1p']['base'] = np.exp(1)    # cosmetic

# Filter too lncRNAs, antisense RNAs, and lowly expressed genes for DE analysis
min_n_cells = adata.obs[chosen].value_counts().min() * .5
test = (~adata.var_names.str.contains('.', regex=False)) & \
       (~adata.var_names.str.contains('LIN')) & \
       (np.sum(adata.layers['raw'].toarray()>0, axis=0) >= min_n_cells)
adata = adata[:,test].copy()

# Marker genes
sc.tl.rank_genes_groups(adata, groupby=chosen, method='wilcoxon', pts=True)

# Filter and isualize marker genes
res_de = mt.ut.format_rank_genes_groups(adata)
order_groups = mt.ut.order_groups(adata, groupby=chosen, obsm_key='X_pca', n_dims=20)
res_de_filtered = mt.ut.format_rank_genes_groups(adata, filter_genes=True)
top_markers = mt.ut.get_top_markers(res_de_filtered, order_groups=order_groups, ntop=3)
df_plot = res_de.query('gene in @top_markers')

fig, ax = plt.subplots(figsize=(6,3.5))
plu.dotplot(df_plot, 'gene', 'group', 
            order_x=top_markers, order_y=order_groups,
            color='log2FC', size='pct_group', ax=ax, vmin=-5, vmax=5)
plu.format_ax(ax=ax, rotx=90, xlabel='', ylabel='Clusters')
ax.get_legend().set_bbox_to_anchor((1,1.2))
fig.tight_layout()
plt.show()

# Annotate individual clusters: e.g. cluster 4 --> IFN response
cluster = "4"
res_de_filtered.query('group==@cluster')
res_gsea = mt.ut.run_GSEA(
    res_de.query('group==@cluster').set_index('gene')['log2FC'], 
    max_size_set=1000, 
    collections='MSigDB_Hallmark_2020'
)[1]
res_gsea['Term'] = res_gsea['Term'].map(lambda x: x.replace('MSigDB_Hallmark_2020__', ''))
res_gsea


# 3. Final cell state annotation --------------------------------------- #

# Final cell state annotation map
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

# Re-perform marker gene selection (after last cell removal)
chosen = 'cell_state'
min_n_cells = adata.obs[chosen].value_counts().min() * .5
test = np.sum(adata.layers['raw'].toarray()>0, axis=0) >= min_n_cells
adata = adata[:,test].copy()
sc.tl.rank_genes_groups(adata, groupby=chosen, method='wilcoxon', pts=True)

# Filter final marker genes
res_de = mt.ut.format_rank_genes_groups(adata)
order_groups = mt.ut.order_groups(adata, groupby=chosen, obsm_key='X_pca', n_dims=20)
res_de_filtered = mt.ut.format_rank_genes_groups(adata, filter_genes=True)
top_markers = mt.ut.get_top_markers(res_de_filtered, order_groups=order_groups, ntop=3)
df_plot = res_de.query('gene in @top_markers')

# Write out filtered markers and pre-processed expression AnnData
res_de_filtered.to_csv(os.path.join(os.path.join(path_results, 'cell_state_markers.csv')))
adata.write(os.path.join(os.path.join(path_data, 'expression.h5ad')))


##
