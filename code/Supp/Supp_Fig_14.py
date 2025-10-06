
"""
Supp Fig 14
--distance_metric parameter tuning. 

# TO DO: Add other non binary distances in the benchmark (e.g. cosine, correlation).
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
path_data = os.path.join(path_main, 'data', 'bench', 'tune_distances')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')


# Set visualization params
plu.set_rcParams({'figure.dpi':350})


##


# Plot metrics with standard or weighted jaccard distance --------------------- # 

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


# Plot
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
fig.savefig(os.path.join(path_figures, 'Supp_Fig_14.pdf'))


##