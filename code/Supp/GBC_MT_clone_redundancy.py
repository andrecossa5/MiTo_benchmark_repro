"""
MT- and GBC- clone redundancy
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
path_expression = os.path.join(path_main, 'data', 'longitudinal')
path_data = os.path.join(path_main, 'results', 'others', 'Fig4', 'longitudinal')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')

# Read data
adata = sc.read(os.path.join(path_expression, 'expression.h5ad'))
afm = sc.read(os.path.join(path_data, 'afm_filtered.h5ad'))
tree_metrics = pd.read_csv(os.path.join(path_data, 'tree_metrics.csv'), index_col=0)
with open(os.path.join(path_data, 'annotated_tree.pickle'), 'rb') as f:
    tree = pickle.load(f)


##


# Get labels and distances
afm.obs['MiTo clone'] = tree.cell_meta['MiTo clone']
df = afm.obs[['MiTo clone', 'GBC']].loc[lambda x: ~x['MiTo clone'].isna()]
D = adata[df.index].obsp['distances'].toarray()

from sklearn.metrics import silhouette_score
silhouette_score(D, labels=df['MiTo clone'])
silhouette_score(D, labels=df['GBC'])

from sklearn.metrics import silhouette_samples
df['sil_MT'] = silhouette_samples(D, labels=df['MiTo clone'])
df['sil_GBC'] = silhouette_samples(D, labels=df['GBC'])

plu.set_rcParams()
fig, ax = plt.subplots(figsize=(3.5, 3.5))
plu.dist(df, x='sil_MT', ax=ax, color='b')
plu.dist(df, x='sil_GBC', ax=ax, color='orange')
plu.format_ax(ax=ax, xlabel='Silhouette score')
plu.add_legend(colors={'MT':'b', 'GBC': 'orange'}, ax=ax, loc='upper left', bbox_to_anchor=(0,1))
fig.tight_layout()
plt.show()

gbc = df['GBC'].value_counts().loc[lambda x: x>=10].index
mt = df['MiTo clone'].value_counts().loc[lambda x: x>=10].index

df_filtered = df.loc[df['MiTo clone'].isin(mt)&df['GBC'].isin(gbc)]
(pd.crosstab(df_filtered['GBC'], df['MiTo clone'])>10).sum(axis=0)


##