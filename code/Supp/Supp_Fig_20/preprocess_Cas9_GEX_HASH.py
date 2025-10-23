"""
Supp Fig 20
- Re-analysis of Weng et. al, 2024 dataset: mouse batch 1.
- Pre-processing of GEX/Cas9 data, after cell demultiplexing cells with hashes library.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import mito as mt
from sklearn.metrics import silhouette_score


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'data', 'redeem_mouse')


##


# De-multiplexing with oligo hashes
hashes = (
    pd.read_csv(
        os.path.join(path_data, 'GSE259284_RAW', 'GSM8113083_Crispr_Mouse_Batch1.Hash_annotation.tsv.gz'), 
        header=None, names=['hash', 'sample'])
    .set_index('hash')
    ['sample'].to_dict()
)
cell_annot = (
    pd.read_csv(os.path.join(path_data, 'GSE259284_RAW', 'GSM8113083_Crispr_Mouse_Batch1.Hash_calls.tsv.gz'), sep='\t')
    .assign(sample=lambda x: x['call'].map(hashes))
    .set_index('cell')
    ['sample']
)


##


# GEX QC
gex = sc.read_10x_h5(
    os.path.join(path_data, 'GSE259284_RAW', 
                 'GSM8113081_Crispr_Mouse_Batch1.filtered_feature_bc_matrix.h5')
)

# Fix cell names
gex.obs_names = gex.obs_names.map(lambda x: x.split('-')[0])
gex.obs_names.isin(cell_annot.index).all()  # Check all cells have hash annotation

# QC
gex.obs['sample'] = cell_annot.loc[gex.obs_names]
gex.obs['nUMIs'] = gex.X.sum(axis=1).A1
gex.obs['detected_genes'] = (gex.X>0).sum(axis=1).A1
gex.obs['MT_perc'] = gex[:,gex.var_names.str.startswith('mt-')].X.sum(axis=1).A1 / gex.obs['nUMIs']

# Cell QC GEX
test = (gex.obs['nUMIs']>=500) & (gex.obs['detected_genes']>=250) & (gex.obs['MT_perc']<=.15)
gex = gex[test,:].copy()


##


# Prep Cas9 allele table
CAS9 = []
files = [ x for x in os.listdir(os.path.join(path_data, 'GSE259284_RAW')) if x.endswith('.allele') ]
for file in files:
    df_ = pd.read_csv(os.path.join(path_data, 'GSE259284_RAW', file), index_col=0)
    df_['sample'] = file.split('.')[0]
    CAS9.append(df_)

allele_table = pd.concat(CAS9)
allele_table = allele_table.dropna()

# Write allele table (only cells present in GEX after QC)
assert allele_table['cellBC'].isin(gex.obs_names).all()
allele_table.to_csv(os.path.join(path_data, 'redeem_mouse_batch1_allele_table.tsv'), sep='\t')


##


# Make each sample cas9 AFMs, and pre-process this modality
samples = ['DKO_1', 'CD47_1', 'Ctrl_2', 'Ctrl_1', 'CD47_2', 'CD24_1']

# Here we go
for sample in samples:

    # Make AFM
    cas9 = mt.io.make_afm(
        path_ch_matrix=os.path.join(path_data, 'redeem_mouse_batch1_allele_table.tsv'), 
        sample=sample,
        scLT_system='Cas9',
        kwargs={'priors_groupby':'sample', 'sample_col':'sample'}
    )

    # Distances, kNN, and embeddings
    mt.pp.compute_distances(cas9, precomputed=True, metric='weighted_hamming')
    idx, dists, conn = mt.pp.kNN_graph(D=cas9.obsp['distances'].toarray(), k=15, from_distances=True)
    mt.pp.reduce_dimensions(cas9, metric='weighted_hamming')

    # Leiden clustering
    LABELS = []
    SILHOUETTES = []
    D = cas9.obsp['distances'].toarray()
    for r in np.linspace(0.5, 2, 10):
        labels = mt.tl.leiden_clustering(conn, res=r)
        sil = silhouette_score(X=D, labels=labels, metric='precomputed')
        LABELS.append(labels)
        SILHOUETTES.append(sil)

    best_r = np.linspace(0.5, 2, 10)[np.argmax(SILHOUETTES)]
    cas9.obs['cas9_clones'] = mt.tl.leiden_clustering(conn, res=best_r)
    cas9.obs['cas9_clones'] = cas9.obs['cas9_clones'].astype('category')

    # Save cas9 AFM
    if not os.path.exists(os.path.join(path_data, sample)):
        os.mkdir(os.path.join(path_data, sample))
    
    # Save
    cas9.write(os.path.join(path_data, sample, 'cas9.h5ad'))


##
