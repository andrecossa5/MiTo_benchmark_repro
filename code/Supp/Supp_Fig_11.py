"""
Supp Fig 11
MiTo genotyping (--bin_method == MiTo) hyper-parameters tuning:
    1. Binomial Mixture probability threshold: --t_prob
    2. Minimum cell prevalence: --min_cell_prevalence
"""

import os
import pandas as pd
import mito as mt
import matplotlib
import matplotlib.pyplot as plt
import plotting_utils as plu
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'data', 'bench', 'tune_genotyping')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')


# Set visualization params
plu.set_rcParams({'figure.dpi':350})


##


# 1. Get metrics  -------------------------- 

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


# 2. Plot ARI across samples, t_prob and min_cell_prevalence parameter values  -------------------------- 

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
fig.savefig(os.path.join(path_figures, 'Supp_Fig_11.pdf'))


##