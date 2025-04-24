"""
Bench distances
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import mito as mt
import matplotlib
import matplotlib.pyplot as plt
import plotting_utils as plu
from itertools import product
from scipy.sparse import csr_matrix
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'data', 'bench', 'tune_distances')
# path_figures = os.path.join(path_main, 'results', 'figures', 'Fig2')
# path_results = os.path.join(path_main, 'results', 'others', 'Fig2')


##


# 1. Extended summary. ---------------------# 

L = []
for folder,_,files in os.walk(path_data):
    if any([ x.startswith('all') for x in files]):
        df,metrics,options = mt.ut.format_tuning(folder)
        L.append(df)

df = pd.concat(L)

varying_options = (df[options].nunique()).loc[lambda x:x>1].index.to_list()
metrics_of_interest = ['ARI', 'NMI', 'corr', 'AUPRC', 'mean_CI', 'mean_RI', 'n MiTo clone']
# (
#     df.groupby(['sample']+varying_options)
#     [['ARI', 'NMI', 'corr', 'AUPRC', 'n_cells', 'n_vars', 'n_GBC_groups']]
#     .median()
# )

##

fig, axs =  plt.subplots(1,len(metrics_of_interest),figsize=(14,5))

x_order = ['MDA_clones', 'MDA_lung', 'MDA_PT']
by_order = ['jaccard', 'weighted_jaccard']

cmap = plu.create_palette(df, 'metric', order=by_order, palette='Reds')
for i,metric in enumerate(metrics_of_interest):
    plu.bar(
        df.query('min_AD=="2"'), 
        x='sample', y=metric, 
        by='metric',
        x_order=x_order,
        by_order=by_order,
        categorical_cmap=cmap,
        ax=axs[i]
    )
    plu.format_ax(ax=axs[i], xlabel='', ylabel=metric, reduced_spines=True, rotx=90)

plu.add_legend(cmap, ax=axs[3], label='Distance metric', loc='center', bbox_to_anchor=(.5,1.2), ncols=2)
fig.subplots_adjust(top=.75, left=.1, bottom=.25, right=.9, wspace=.6)
plt.show()


##

## 2. Perturbation experiment, and corr/ARI distances.  ---------------------# 

path_data = os.path.join(path_main, 'data', 'general', 'AFMs', 'maegatk')

# Read and filter afm and cells
sample = 'MDA_clones'
afm = sc.read(os.path.join(path_data, f'afm_{sample}.h5ad'))
afm = mt.pp.filter_cells(afm, cell_filter='filter2')
afm = mt.pp.filter_afm(afm, lineage_column='GBC', min_cell_number=5)

mt.pp.compute_distances(afm, precomputed=True, distance_key='jaccard', metric='jaccard')
tree_jacc = mt.tl.build_tree(afm, distance_key='jaccard', precomputed=True)
corr_jacc = mt.ut.calculate_corr_distances(tree_jacc)[0]

mt.pp.compute_distances(afm, precomputed=True, distance_key='weighted_jaccard', metric='weighted_jaccard')
tree_wjacc = mt.tl.build_tree(afm, distance_key='weighted_jaccard', precomputed=True)
corr_wjacc = mt.ut.calculate_corr_distances(tree_wjacc)[0]

fig, axs = plt.subplots(1,2,figsize=(8,4))
mt.pl.heatmap_distances(afm, distance_key='jaccard', tree=tree_jacc, ax=axs[0])
mt.pl.heatmap_distances(afm, distance_key='weighted_jaccard', tree=tree_wjacc, ax=axs[1])
fig.tight_layout()
plt.show()


##


L = []
combos = product(np.linspace(.1,1,5), np.linspace(.1,50,5), [True, False])

for i,(perc_sites,theta,add) in enumerate(combos):

    stats = {}
    stats['perc_sites'] = perc_sites
    stats['theta'] = theta
    stats['add'] = add
    afm_perturbed, AD_corr = mt.ut.perturb_AD_counts(afm.copy(), perc_sites=perc_sites, theta=theta, add=add)
    stats['AD_corr'] = AD_corr

    for metric in ['jaccard', 'weighted_jaccard']:
        d_ = stats.copy()
        d_['metric'] = metric
        mt.pp.call_genotypes(afm_perturbed, bin_method='MiTo')
        mt.pp.compute_distances(afm_perturbed, precomputed=True, distance_key=metric, metric=metric)
        index,_,_ = mt.pp.kNN_graph(D=afm_perturbed.obsp[metric].A, from_distances=True, k=10)
        d_['NN_purity'] = mt.ut.NN_purity(index, labels=afm_perturbed.obs['GBC'].values)
        d_['kBET'] = mt.ut.kbet(index, batch=afm_perturbed.obs['GBC'])
        tree_perturbed = mt.tl.build_tree(afm_perturbed, distance_key=metric, precomputed=True)
        d_['tree_char_corr'] = mt.ut.calculate_corr_distances(tree_perturbed)[0]
        L.append(d_)


##

df_plot = pd.DataFrame(L)
df_plot['perc_sites'] = np.round(df_plot['perc_sites']*100,1)
df_plot['theta'] = np.round(df_plot['theta'],1)

# Viz
fig, axs = plt.subplots(2,2,figsize=(8,5))


# Viz
by_order = ['jaccard', 'weighted_jaccard']
metrics = ['kBET', 'tree_char_corr']
names = dict(zip(metrics, ['kBET', 'Tree-character\ndistance correlation']))
cmap = plu.create_palette(df_plot, 'metric', order=by_order, palette='Reds')

for i,metric in enumerate(metrics):
    plu.strip(df_plot, x='perc_sites', y=metric, ax=axs[0,i], by='metric', categorical_cmap=cmap, by_order=by_order, size=2)
    plu.box(df_plot, x='perc_sites', y=metric, ax=axs[0,i], by='metric', categorical_cmap=cmap, by_order=by_order)
    plu.format_ax(ax=axs[0,i], reduced_spines=True, xlabel='% perturbed sites', ylabel=names[metric])
    
for i,metric in enumerate(metrics):
    plu.strip(df_plot, x='theta', y=metric, ax=axs[1,i], by='metric', categorical_cmap=cmap, by_order=by_order, size=2)
    plu.box(df_plot, x='theta', y=metric, ax=axs[1,i], by='metric', categorical_cmap=cmap, by_order=by_order)
    plu.format_ax(ax=axs[1,i], reduced_spines=True, xlabel='Binomial noise', ylabel=names[metric])

fig.tight_layout()
plt.show()





