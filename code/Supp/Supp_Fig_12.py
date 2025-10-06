"""
Supp Fig 12
MiTo genotyping --bin_method benchmark vs "vanilla" baseline method.
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


# 1. Get metrics  -------------------------- #

# Format metrics df 
L = []
for folder,_,files in os.walk(path_data):
    if any([ x.startswith('all') for x in files]):
        df,metrics,options = mt.ut.format_tuning(folder)
        L.append(df)
df = pd.concat(L)
df.loc[df['pp_method'] == 'mito_preprocessing', 'pp_method'] = 'MiTo'

varying_options = (df[options].nunique()).loc[lambda x:x>1].index.to_list()
metrics_of_interest = ['ARI', 'NMI', 'corr', 'AUPRC', 'n_cells', 'n_vars', 'n_GBC_groups']
metrics_of_interest += [ 'median_n_vars_per_cell' ]


##


# 2. Metrics of interest: MiTo vs vanilla genotyping -------------------------- #

cmap = {'vanilla' : '#E8E0E0', 'MiTo' : '#D15757' }

for sample in ['MDA_clones', 'MDA_PT', 'MDA_lung']:

    # Fig
    fig, axs = plt.subplots(1,len(metrics_of_interest), figsize=(12,2.5))

    for i,metric in enumerate(metrics_of_interest):

        if sample == 'MDA_PT':
            df_mito = df.query('sample==@sample and bin_method=="MiTo" and t_prob=="0.7" and min_cell_prevalence=="0.05"')  
        else:
            df_mito = df.query('sample=="MDA_PT" and bin_method=="MiTo" and t_prob=="0.7" and min_cell_prevalence=="0.1"')

        df_vanilla = df.query('sample==@sample and bin_method=="vanilla"')
        df = pd.concat([df_vanilla,df_mito]).set_index('job_id')

        plu.bar(
            df,
            x='min_AD', y=metric, by='bin_method', ax=axs[i], 
            categorical_cmap=cmap, by_order=['vanilla', 'MiTo'], x_order=["1","2"]
        )
        plu.format_ax(ax=axs[i], reduced_spines=True, xlabel='min n alt UMIs', ylabel=metric)

    fig.tight_layout()
    fig.savefig(os.path.join(path_figures, f'Supp_Fig_12_{sample}.pdf'))


##