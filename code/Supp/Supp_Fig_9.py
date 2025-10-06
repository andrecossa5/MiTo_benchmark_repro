"""
Supp Fig 9
Tuning of MiTo MT-SNVs filter hyper-parameters (see nf-MiTo Allele Frequency Matrix Preprocessing parameters)
"""

import os
import pandas as pd
import mito as mt
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import plotting_utils as plu
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_afm = os.path.join(path_main, 'data', 'general', 'AFMs')
path_tuning = os.path.join(path_main, 'data', 'bench', 'tune_filtering')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')


# Set visualization params
plu.set_rcParams({'figure.dpi':350})


##


# 1 Read data -----------------------------------#

df_all_filter,_,_ = mt.ut.format_tuning(path_tuning)
df_all_filter = df_all_filter.assign(optional='DBs and spatial')
df_no_dbs,_,_ = mt.ut.format_tuning(os.path.join(path_tuning, 'no_dbs'))
df_no_dbs = df_no_dbs.assign(optional='no DBs')
df_no_spatial,metrics,options = mt.ut.format_tuning(os.path.join(path_tuning, 'no_spatial'))
df_no_spatial = df_no_spatial.assign(optional='no spatial')
df = pd.concat([df_all_filter, df_no_dbs, df_no_spatial])

# Format options
varying_options = ['sample', 'af_confident_detection', 'min_n_confidently_detected', 'min_mean_AD_in_positives', 'optional']
df = df.rename(columns={'average_degree':'Av. degree', 'transitions_vs_transversions_ratio':'Trasitions /\ntransversions'})
metrics_of_interest = ['ARI', 'n_cells', 'n_vars', 'Av. degree', 'Trasitions /\ntransversions']


# 2. Supp Fig 9a: Allelic Frequency spectrum of unfiltered AFMs -----------------------------------#

fig, axs = plt.subplots(1,3,figsize=(7.5,3), sharey=True)

samples = ['MDA_clones', 'MDA_PT', 'MDA_lung']
for i,sample in enumerate(samples):
    afm = sc.read(os.path.join(path_afm, sample, 'afm_unfiltered.h5ad'))
    afm = mt.pp.filter_cells(afm, cell_filter='filter2')
    mt.pl.vars_AF_spectrum(afm, ax=axs[i], color='k', linewidth=.1)
    axs[i].set(title=sample, ylabel='Allelic Frequency' if i==0 else '')

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'Supp_Fig_9a.pdf'))


# 3. Supp Fig 9b: n of MT-SNVs excluded by matching dbSNPs or REDIdb databases queries -----------------------------------#

fig, axs = plt.subplots(2,1,figsize=(2.5,3), sharex=True)

order = ['MDA_clones', 'MDA_PT', 'MDA_lung']
plu.strip(df.set_index('job_id'), x='sample', y='n_dbSNP', ax=axs[0], x_order=order, linewidth=.5, color='white')
plu.format_ax(ax=axs[0], xlabel='', ylabel='n dbSNP MT-SNVs', reduced_spines=True)
plu.strip(df.set_index('job_id'), x='sample', y='n_REDIdb', ax=axs[1], x_order=order, linewidth=.5, color='white')
plu.format_ax(ax=axs[1], xlabel='', ylabel='n REDIdb MT-SNVs', reduced_spines=True)

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'Supp_Fig_9b.pdf'))


# 2. Supp Fig 9c-f: Extended summary -----------------------------------#

parameters = [
    'min_mean_AD_in_positives', 
    'af_confident_detection', 
    'min_n_confidently_detected', 
    'optional'
]
subtitles = [
    'min AD in +cells', 
    'AF confident detection', 
    'Min confident cells', 
    'Optional filters'
]
subfigure_idxs = ['c','d','e','f']
palettes = ['Reds', 'Blues', 'Greens', 'Oranges']


for parameter,subtitle,subfigure_idx,palette in zip(parameters,subtitles,subfigure_idxs,palettes):

    fig, axs = plt.subplots(5,1,figsize=(4,9), sharex=True)

    x_order = ['MDA_clones', 'MDA_PT', 'MDA_lung']

    # Handle parameter
    if parameter == 'optional':
        by_order = ['no DBs', 'no spatial', 'DBs and spatial']
    else:
        try:
            to_numeric = df[parameter].astype('int')
        except:
            to_numeric = df[parameter].astype('float')
        by_order = [ str(x) for x in sorted(to_numeric.unique()) ]
    
    # Viz
    cmap = plu.create_palette(df, parameter, order=by_order, palette=palette)
    df_ = df.set_index('job_id') # Reindex for plotting

    for i,metric in enumerate(metrics_of_interest):
        plu.bar(
            df_, 
            x='sample', 
            y=metric, 
            x_order=x_order, 
            by=parameter, 
            by_order=by_order, 
            ax=axs[i], 
            categorical_cmap=cmap
        )
        plu.format_ax(ax=axs[i], ylabel=metric, xlabel='', reduced_spines=True)

    plu.add_legend(cmap, label=subtitle, ax=axs[0], bbox_to_anchor=(0.5, 1.3), loc='center', ncols=3)
    fig.subplots_adjust(left=.2, right=.85, bottom=.1, top=.9)
    fig.savefig(os.path.join(path_figures, f'Supp_Fig_9{subfigure_idx}.pdf'))


##