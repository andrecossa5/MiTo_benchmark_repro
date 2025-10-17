"""
Supp Fig 15
- Re-analysis of Weng et. al, 2024 dataset: mouse batch 1.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import cassiopeia as cs
import mito as mt
import matplotlib
import matplotlib.pyplot as plt
import plotting_utils as plu
from scipy.sparse import csr_matrix
from anndata import AnnData
matplotlib.use('macOSX')


##


# Set paths
# path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
# path_data = os.path.join(path_main, 'data', 'ludwig_2019')
# path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')
path_data = '/Users/IEO5505/Desktop/MI_TO/last_data/redeem_mouse/GSE259284_RAW'

# Set visualization params
plu.set_rcParams({'figure.dpi':100})


##


# De-multiplexing with oligo hashes
hashes = (
    pd.read_csv(
        os.path.join(path_data, 'GSM8113083_Crispr_Mouse_Batch1.Hash_annotation.tsv.gz'), 
        header=None, names=['hash', 'sample'])
    .set_index('hash')
    ['sample'].to_dict()
)
cell_annot = (
    pd.read_csv(os.path.join(path_data, 'GSM8113083_Crispr_Mouse_Batch1.Hash_calls.tsv.gz'), sep='\t')
    .assign(sample=lambda x: x['call'].map(hashes))
    .set_index('cell')
    ['sample']
)
cell_annot.index = cell_annot.index.map(lambda x: f'{x}-1')

# GEX
gex = sc.read_10x_h5(os.path.join(path_data, "GSM8113081_Crispr_Mouse_Batch1.filtered_feature_bc_matrix.h5"))
gex.var
gex.obs_names.isin(cell_annot.index).all()
gex.obs['sample'] = cell_annot.loc[gex.obs_names]
gex.obs['nUMIs'] = gex.X.sum(axis=1).A1
gex.obs['detected_genes'] = (gex.X>0).sum(axis=1).A1
gex.obs['MT_perc'] = gex[:,gex.var_names.str.startswith('mt-')].X.sum(axis=1).A1 / gex.obs['nUMIs']

# Cell QC GEX
test = (gex.obs['nUMIs']>=500) & (gex.obs['detected_genes']>=250) & (gex.obs['MT_perc']<=.15)
gex = gex[test,:].copy()


##


# Pre-process Cas9

# Load allele tables
CAS9 = []
files = [ x for x in os.listdir(path_data) if x.endswith('.allele') ]
for file in files:
    df_ = pd.read_csv(os.path.join(path_data, file), index_col=0)
    df_['sample'] = file.split('.')[0]
    CAS9.append(df_)

# Priors 
allele_table = pd.concat(CAS9)
allele_table = allele_table.dropna()
indel_priors = cs.pp.compute_empirical_indel_priors(
    allele_table, grouping_variables=["intBC", "sample"]
)

# Character matrix
(
    char_matrix,
    priors,
    _,
) = cs.pp.convert_alleletable_to_character_matrix(
    allele_table, allele_rep_thresh=0.9, mutation_priors=indel_priors
)
char_matrix.index = char_matrix.index.map(lambda x: f'{x}-1')

# Subset GEX and Cas9 char matrix
cells = list(set(gex.obs_names) & set(char_matrix.index))
gex = gex[cells,:].copy()
char_matrix = char_matrix.loc[cells]

# Cas9 AFM (MiTo format)
gex.obs['sample'].value_counts()
sample = 'CD47_1' # One sample
cells = gex.obs.query('sample==@sample').index

# AFM
from mito.io.format_afm import _add_priors
cas9 = AnnData(
    X=csr_matrix(char_matrix.loc[cells].values), 
    obs=gex.obs.loc[cells], 
    var=pd.DataFrame(index=char_matrix.columns),
    uns={
       'pp_method':'cassiopeia', 
       'scLT_system':'Cas9', 
    }
)
_add_priors(cas9, priors)
cas9.layers['bin'] = cas9.X.copy()

# Distances, kNN, and embeddings
mt.pp.compute_distances(cas9, precomputed=True, metric='weighted_hamming')
idx, dists, conn = mt.pp.kNN_graph(D=cas9.obsp['distances'].toarray(), k=15, from_distances=True)
mt.pp.reduce_dimensions(cas9, metric='weighted_hamming')

# Leiden clustering
from sklearn.metrics import silhouette_score

LABELS = []
SILHOUETTES = []
D = cas9.obsp['distances'].toarray()
for r in np.linspace(0.5, 2, 10):
    labels = mt.tl.leiden_clustering(conn, res=r)
    sil = silhouette_score(X=D, labels=labels, metric='precomputed')
    LABELS.append(labels)
    SILHOUETTES.append(sil)

plt.plot(np.linspace(0.5, 2, 10), SILHOUETTES)
plt.show()

best_r = np.linspace(0.5, 2, 10)[np.argmax(SILHOUETTES)]
cas9.obs['cas9_clones'] = mt.tl.leiden_clustering(conn, res=best_r)
cas9.obs['cas9_clones'] = cas9.obs['cas9_clones'].astype('category')

fig, ax = plt.subplots(figsize=(3,3))
sc.pl.embedding(cas9, basis='X_umap', color='cas9_clones', ax=ax, show=False, frameon=False, size=75)
fig.tight_layout()
plt.show()


##


# Assemble MiTo AFM
basecalls = pd.read_csv(os.path.join(path_data, 'filtered_basecalls.csv'), index_col=0)
basecalls.columns = ['genotype', 'cell', 'mut', 'AD', 'DP', 'type', 'context', 'AD_before_trim', 'AF']
basecalls['pos'] = basecalls['mut'].map(lambda x: x.split('_')[0])
basecalls['type'] = basecalls['type'].map(lambda x: x.replace('_', '>'))
basecalls = basecalls.drop(columns=['genotype']).assign(mut=lambda x: x['pos']+'_'+x['type'])
basecalls = basecalls[['cell', 'mut', 'AD', 'DP', 'AF']]

## Raw cell-MT_genome position consensus UMI counts
cov = pd.read_csv(os.path.join(path_data, 'coverage.csv'))
cov.columns = ['cell', 'pos', 'total', 'less stringent', 'stringent', 'very stringent']
cov = (
    ## Retain total coverage consensus UMI counts for DP values,
    ## as in https://github.com/chenweng1991/redeemR/blob/master/R/VariantSummary.R
    cov[['cell', 'pos', 'total']].rename(columns={'total':'DP'})        
)

# Pivot filtered basecalls, and filter cell_meta and cov for cells
AD = basecalls.pivot(index='cell', columns='mut', values='AD').fillna(0).copy()
DP = basecalls.pivot(index='cell', columns='mut', values='DP').fillna(0).copy()

# Checks
assert (AD.index.value_counts()==1).all()
assert (np.sum(DP>0, axis=1)>0).all()
assert (np.sum(DP>0, axis=1)>0).all() 
assert (AD.index == DP.index).all()
assert (AD.columns == DP.columns).all()

# Char and cell metadata
char_meta = DP.columns.to_series().to_frame('mut')
char_meta['pos'] = char_meta['mut'].map(lambda x: int(x.split('_')[0]))
char_meta['ref'] = char_meta['mut'].map(lambda x: x.split('_')[1].split('>')[0])
char_meta['alt'] = char_meta['mut'].map(lambda x: x.split('_')[1].split('>')[1])
char_meta = char_meta[['pos', 'ref', 'alt']]

# Create sparse matrices, and store into the AnnData object
AF = csr_matrix(np.divide(AD.values,(DP.values+.00000001)).astype(np.float32))
AD = csr_matrix(AD.values).astype(np.int16)
DP = csr_matrix(DP.values).astype(np.int16)
afm = AnnData(
    X=AF, 
    obs=cell_meta, 
    var=char_meta, 
    layers={'AD':AD, 'DP':DP}, 
    uns={'pp_method':pp_method, 'scLT_system':scLT_system, 'raw_basecalls_metrics':None}
)
sorted_vars = afm.var['pos'].sort_values().index
assert sorted_vars.size == afm.shape[1]
afm = afm[:,sorted_vars].copy()

# Add complete site coverage info
cov = cov.query('cell in @filtered_cells').pivot(index='cell', columns='pos', values='DP').fillna(0)
mapping = afm.var['pos'].to_dict()
df_ = pd.DataFrame({ mut : cov[mapping[mut]].values for mut in mapping }, index=filtered_cells)
assert all(df_.columns == afm.var_names)
afm.layers['site_coverage'] = csr_matrix(df_.values)
afm.obs['mean_site_coverage'] = cov.mean(axis=1)   