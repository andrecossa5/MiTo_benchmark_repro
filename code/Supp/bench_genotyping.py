"""
Bench genotyping
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import mito as mt
import matplotlib
import matplotlib.pyplot as plt
import plotting_utils as plu
from scipy.sparse import csr_matrix
from itertools import product
matplotlib.use('macOSX')



##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'data', 'bench', 'tune_genotyping')
# path_figures = os.path.join(path_main, 'results', 'figures', 'Fig2')
# path_results = os.path.join(path_main, 'results', 'others', 'Fig2')


##


# Params
plu.set_rcParams()
matplotlib.rcParams.update({'figure.dpi':150})


##


# 1. Tune MiTo hyperparameters -------------------------- 

# Format metrics df 
L = []
for folder,_,files in os.walk(path_data):
    if any([ x.startswith('all') for x in files]):
        df,metrics,options = mt.ut.format_tuning(folder)
        L.append(df)
df = pd.concat(L)
df.loc[df['pp_method'] == 'mito_preprocessing', 'pp_method'] = 'MiTo'

varying_options = (df[options].nunique()).loc[lambda x:x>1].index.to_list()
metrics_of_interest = ['ARI', 'corr', 'n_cells', 'n_GBC_groups']


##

fig, axs = plt.subplots(1,3,figsize=(12,2.5))

by_order = ['0.01', '0.05', '0.1']
x_order = ['0.5', '0.7', '0.9']
cmap = plu.create_palette(df, 'min_cell_prevalence', palette='Oranges', order=by_order)

for i,sample in enumerate(['MDA_clones', 'MDA_PT', 'MDA_lung']):
    plu.bar(
        df.query('bin_method=="MiTo" and sample==@sample'), 
        ax=axs[i], x='t_prob', y='ARI', by='min_cell_prevalence',
        x_order=x_order, by_order=by_order, categorical_cmap=cmap
    )
    plu.format_ax(ax=axs[i], ylabel='ARI', title=sample,
        xlabel='Binomial mixture \nprobability threshold', reduced_spines=True)

plu.add_legend(cmap, label='Min prevalence', ax=axs[i])
fig.subplots_adjust(left=.1, right=.8, top=.9, bottom=.25)
plt.show()


##


# 2. vanilla vs top MiTo, one per sample -------------------------- 

# Define options and metrics
metrics_of_interest = ['ARI', 'NMI', 'corr', 'AUPRC', 'n_cells', 'n_vars', 'n_GBC_groups']
metrics_of_interest += [ 'median_n_vars_per_cell' ]

# Separate sample image
fig, axs = plt.subplots(1,len(metrics_of_interest), figsize=(12,2.5))

sample = 'MDA_PT'
cmap = {'vanilla' : '#E8E0E0', 'MiTo' : '#D15757' }

for i,metric in enumerate(metrics_of_interest):
    plu.bar(
        pd.concat(
            [
                df.query('sample==@sample and bin_method=="vanilla"'),
                df.query('sample==@sample and bin_method=="MiTo" and t_prob=="0.7" and min_cell_prevalence=="0.1"'),
                # df.query('sample=="MDA_PT" and bin_method=="MiTo" and t_prob=="0.7" and min_cell_prevalence=="0.05"')
            ]
        ).set_index('job_id'), 
        x='min_AD', y=metric, by='bin_method', ax=axs[i], 
        categorical_cmap=cmap, by_order=['vanilla', 'MiTo'], x_order=["1","2"]
    )
    plu.format_ax(ax=axs[i], reduced_spines=True, xlabel='min n alt UMIs', ylabel=metric)

fig.tight_layout()
plt.show()


##


# 3. Drop-out and spurious counts experiment -------------------------- 

path_data = os.path.join(path_main, 'data', 'general', 'AFMs', 'maegatk')

# Read and filter afm and cells
sample = 'MDA_clones'
afm = sc.read(os.path.join(path_data, f'afm_{sample}.h5ad'))
afm = mt.pp.filter_cells(afm, cell_filter='filter2')
afm = mt.pp.filter_afm(afm, lineage_column='GBC', min_cell_number=5)


##


L = []
combos = product(np.linspace(.1,1,10), np.linspace(.1,50,10), [True, False])

for i,(perc_sites,theta,add) in enumerate(combos):

    stats = {}
    stats['perc_sites'] = perc_sites
    stats['theta'] = theta
    stats['add'] = add
    afm_perturbed, AD_corr = mt.ut.perturb_AD_counts(afm.copy(), perc_sites=perc_sites, theta=theta, add=add)
    stats['AD_corr'] = AD_corr

    for bin_method in ['MiTo', 'vanilla']:
        d_ = stats.copy()
        d_['bin_method'] = bin_method
        mt.pp.call_genotypes(afm_perturbed, bin_method=bin_method)
        mt.pp.compute_distances(afm_perturbed, precomputed=True, distance_key='perturbed')
        index,_,_ = mt.pp.kNN_graph(D=afm_perturbed.obsp['perturbed'].A, from_distances=True, k=10)
        d_['NN_purity'] = mt.ut.NN_purity(index, labels=afm_perturbed.obs['GBC'].values)
        d_['kBET'] = mt.ut.kbet(index, batch=afm_perturbed.obs['GBC'])
        tree_perturbed = mt.tl.build_tree(afm_perturbed, distance_key='perturbed', precomputed=True)
        d_['tree_char_corr'] = mt.ut.calculate_corr_distances(tree_perturbed)[0]
        L.append(d_)

##

df = pd.DataFrame(L)
df['perc_sites'] = np.round(df['perc_sites']*100,1)
df['theta'] = np.round(df['theta'],1)

# Viz
fig, axs = plt.subplots(2,2,figsize=(8,5))

by_order = ['vanilla', 'MiTo']
metrics = ['kBET', 'tree_char_corr']
names = dict(zip(metrics, ['kBET', 'Tree-character\ndistance correlation']))
cmap = {'vanilla' : '#E8E0E0', 'MiTo' : '#D15757' }

for i,metric in enumerate(metrics):
    plu.strip(df, x='perc_sites', y=metric, ax=axs[0,i], by='bin_method', categorical_cmap=cmap, by_order=by_order, size=2)
    plu.box(df, x='perc_sites', y=metric, ax=axs[0,i], by='bin_method', categorical_cmap=cmap, by_order=by_order)
    plu.format_ax(ax=axs[0,i], reduced_spines=True, xlabel='% perturbed sites', ylabel=names[metric])
    
for i,metric in enumerate(metrics):
    plu.strip(df, x='theta', y=metric, ax=axs[1,i], by='bin_method', categorical_cmap=cmap, by_order=by_order, size=2)
    plu.box(df, x='theta', y=metric, ax=axs[1,i], by='bin_method', categorical_cmap=cmap, by_order=by_order)
    plu.format_ax(ax=axs[1,i], reduced_spines=True, xlabel='Binomial noise', ylabel=names[metric])

fig.tight_layout()
plt.show()


##


fig, ax = plt.subplots()
plu.add_legend(cmap, label='Genotyping', ax=ax, loc='center', bbox_to_anchor=(.5,.5), ncols=2)
fig.tight_layout()
plt.show()


##