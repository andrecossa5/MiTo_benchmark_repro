
"""
Supp Fig 13
Simulation experiment. Drop-out and error noise.
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
matplotlib.use('macOSX')



##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'data', 'general', 'AFMs')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')


# Set visualization params
plu.set_rcParams({'figure.dpi':350})


# 1. Simulate genotyping with different levels of drop-out and sequencing noise -------------------------- #

# Read and filter AFM
sample = 'MDA_clones'
afm = sc.read(os.path.join(path_data, sample, f'afm_unfiltered.h5ad'))
afm = mt.pp.filter_cells(afm, cell_filter='filter2')
afm = mt.pp.filter_afm(afm, lineage_column='GBC', min_cell_number=1)


##


L = []
combos = list(product(np.linspace(.1,1,5), np.linspace(.01,5,5), [True]))

for step in range(50):

    for i,(perc_sites,theta,add) in enumerate(combos):

        stats = {}
        stats['step'] = step
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
            d_['AUPRC'] =  mt.ut.distance_AUPRC(afm_perturbed.obsp['perturbed'].toarray(), afm_perturbed.obs['GBC'].values)
            tree_perturbed = mt.tl.build_tree(afm_perturbed, distance_key='perturbed', precomputed=True)
            d_['tree_char_corr'] = mt.ut.calculate_corr_distances(tree_perturbed)[0]
            L.append(d_)

##


# 2. Visualize simulation results -------------------------- #

df = pd.DataFrame(L)
df['perc_sites'] = np.round(df['perc_sites']*100)
df['theta'] = np.round(df['theta'],2)
plt.plot(df['perc_sites'], df['AD_corr'], 'ko')
plt.show()


##


# Plot
fig, axs = plt.subplots(2,2,figsize=(8,5))

by_order = ['vanilla', 'MiTo']
metrics = ['AUPRC', 'tree_char_corr']
names = dict(zip(metrics, ['AUPRC', 'Tree-character\ndistance correlation']))
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
fig.savefig(os.path.join(path_figures, 'Supp_Fig_13.pdf'))


##