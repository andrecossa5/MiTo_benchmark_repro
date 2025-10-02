"""
Bench filtering.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import mito as mt
import matplotlib
import matplotlib.pyplot as plt
import plotting_utils as plu
from mito.pp.filters import mask_mt_sites
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'data', 'bench', 'tune_filtering')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')
# path_results = os.path.join(path_main, 'results', 'others', 'Fig2')


##


# Read data
df_all_filter,_,_ = mt.ut.format_tuning(path_data)
df_all_filter = df_all_filter.assign(optional='DBs and spatial')
df_no_dbs,_,_ = mt.ut.format_tuning(os.path.join(path_data, 'no_dbs'))
df_no_dbs = df_no_dbs.assign(optional='no DBs')
df_no_spatial,metrics,options = mt.ut.format_tuning(os.path.join(path_data, 'no_spatial'))
df_no_spatial = df_no_spatial.assign(optional='no spatial')
df = pd.concat([df_all_filter, df_no_dbs, df_no_spatial])

# Format options
varying_options = ['sample', 'af_confident_detection', 'min_n_confidently_detected', 'min_mean_AD_in_positives', 'optional']
df = df.rename(columns={'average_degree':'Av. degree', 'transitions_vs_transversions_ratio':'Trasitions /\ntransversions'})
metrics_of_interest = ['ARI', 'n_cells', 'n_vars', 'Av. degree', 'Trasitions /\ntransversions']


##


# Viz
plu.set_rcParams()
matplotlib.rcParams.update({'figure.dpi':350})

##


# 1. Extended summary -----------------------------------#

fig, axs = plt.subplots(5,1,figsize=(4,9), sharex=True)

x_order = ['MDA_clones', 'MDA_PT', 'MDA_lung']
by_order = [ str(x) for x in np.sort(df['min_mean_AD_in_positives'].astype(float).unique()) ]
cmap = plu.create_palette(df, 'min_mean_AD_in_positives', order=by_order, palette='Reds')

for i,metric in enumerate(metrics_of_interest):
    plu.bar(df.set_index('job_id'), x='sample', y=metric, x_order=x_order, 
              by='min_mean_AD_in_positives', by_order=by_order, ax=axs[i], categorical_cmap=cmap)
    plu.format_ax(ax=axs[i], ylabel=metric, xlabel='', reduced_spines=True)

plu.add_legend(cmap, label='min AD in +cells', ax=axs[0], bbox_to_anchor=(0.5, 1.3), loc='center', ncols=3)
fig.subplots_adjust(left=.2, right=.85, bottom=.1, top=.9)
plt.show()

##

fig, axs = plt.subplots(5,1,figsize=(4,9), sharex=True)

x_order = ['MDA_clones', 'MDA_PT', 'MDA_lung']
by_order = [ str(x) for x in np.sort(df['af_confident_detection'].astype(float).unique()) ]
cmap = plu.create_palette(df, 'af_confident_detection', order=by_order, palette='Blues')

for i,metric in enumerate(metrics_of_interest):
    plu.bar(df.set_index('job_id'), x='sample', y=metric, x_order=x_order, 
              by='af_confident_detection', by_order=by_order, ax=axs[i], categorical_cmap=cmap)
    plu.format_ax(ax=axs[i], ylabel=metric, xlabel='', reduced_spines=True)

plu.add_legend(cmap, label='AF confident detection', ax=axs[0], bbox_to_anchor=(0.5, 1.3), loc='center', ncols=2)

fig.subplots_adjust(left=.2, right=.85, bottom=.1, top=.9)
plt.show()

##

fig, axs = plt.subplots(5,1,figsize=(4,9), sharex=True)

x_order = ['MDA_clones', 'MDA_PT', 'MDA_lung']
by_order = [ str(x) for x in np.sort(df['min_n_confidently_detected'].astype(int).unique()) ]
cmap = plu.create_palette(df, 'min_n_confidently_detected', order=by_order, palette='Greens')

for i,metric in enumerate(metrics_of_interest):
    plu.bar(df.set_index('job_id'), x='sample', y=metric, x_order=x_order, 
              by='min_n_confidently_detected', by_order=by_order, ax=axs[i], categorical_cmap=cmap)
    plu.format_ax(ax=axs[i], ylabel=metric, xlabel='', reduced_spines=True)

plu.add_legend(cmap, label='min n confident cells', ax=axs[0], bbox_to_anchor=(0.5, 1.3), loc='center', ncols=2)

fig.subplots_adjust(left=.2, right=.85, bottom=.1, top=.9)
plt.show()

##

fig, axs = plt.subplots(5,1,figsize=(4,9), sharex=True)

x_order = ['MDA_clones', 'MDA_PT', 'MDA_lung']
by_order = ['no DBs', 'no spatial', 'DBs and spatial']
cmap = plu.create_palette(df, 'optional', order=by_order, palette='Oranges')

for i,metric in enumerate(metrics_of_interest):
    plu.bar(df.set_index('job_id'), x='sample', y=metric, x_order=x_order, 
              by='optional', by_order=by_order, ax=axs[i], categorical_cmap=cmap)
    plu.format_ax(ax=axs[i], ylabel=metric, xlabel='', reduced_spines=True)

plu.add_legend(cmap, label='Optional filters', ax=axs[0], bbox_to_anchor=(0.5, 1.3), loc='center', ncols=2)

fig.subplots_adjust(left=.2, right=.85, bottom=.1, top=.9)
plt.show()


##


# 2. n dbs excluded MUTs -----------------------------------#

fig, axs = plt.subplots(2,1,figsize=(2.5,3), sharex=True)

order = ['MDA_clones', 'MDA_PT', 'MDA_lung']
plu.strip(df.set_index('job_id'), x='sample', y='n_dbSNP', ax=axs[0], x_order=order, linewidth=.5, color='white')
plu.format_ax(ax=axs[0], xlabel='', ylabel='n dbSNP MT-SNVs', reduced_spines=True)
plu.strip(df.set_index('job_id'), x='sample', y='n_REDIdb', ax=axs[1], x_order=order, linewidth=.5, color='white')
plu.format_ax(ax=axs[1], xlabel='', ylabel='n REDIdb MT-SNVs', reduced_spines=True)

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'Supp_Fig_9_uninformative.pdf'))


##


# 3. VG like spectrum -----------------------------------#

path_data = os.path.join(path_main, 'data', 'general', 'AFMs', 'maegatk')

fig, axs = plt.subplots(1,3,figsize=(7.5,3), sharey=True)

samples = ['MDA_clones', 'MDA_PT', 'MDA_lung']
for i,sample in enumerate(order):
    afm = sc.read(os.path.join(path_data, f'afm_{sample}.h5ad'))
    afm = mt.pp.filter_cells(afm, cell_filter='filter2')
    mt.pl.vars_AF_spectrum(afm, ax=axs[i], color='k', linewidth=.1)
    axs[i].set(title=sample, ylabel='Allelic Frequency' if i==0 else '')

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'Supp_Fig_9_spectrum.pdf'))


##


# 4. MT-SNVs per clone -----------------------------------#

L = []
fdr_treshold = .05
samples = ['MDA_clones', 'MDA_PT', 'MDA_lung']
for sample in samples:
    afm = sc.read(os.path.join(path_data, f'afm_{sample}.h5ad'))
    afm = mt.pp.filter_cells(afm, cell_filter='filter2')
    afm = mt.pp.filter_afm(afm, lineage_column='GBC', min_cell_number=5, compute_enrichment=True)
    fdrs = afm.var.loc[:,afm.var.columns.str.startswith('FDR')]
    df = pd.DataFrame({
        'cat' : np.where(np.sum(fdrs<=fdr_treshold)>0, 'Defined', 'Undefined'),
        'sample' : [ sample for _ in range(fdrs.columns.size) ],
        'n_muts' : np.sum(fdrs<=fdr_treshold)
    })
    L.append(df)

df_plot = pd.concat(L).reset_index(drop=True)
df_plot['sample'] = pd.Categorical(df_plot['sample'], categories=samples[::-1])

##

fig, ax = plt.subplots(figsize=(4,3))

cmap = plu.create_palette(df_plot, 'cat', order=['Undefined', 'Defined'], palette='Greys')
plu.bb_plot(df_plot, cov1='sample', cov2='cat', ax=ax, categorical_cmap=cmap, show_y=True)
plu.add_legend(cmap, label='Clone type', loc='center', bbox_to_anchor=(.5,1.4), ncols=2, ax=ax)
fig.subplots_adjust(top=.6, bottom=.3, right=.85, left=.25)
fig.savefig(os.path.join(path_figures, 'Supp_Fig_10a.pdf'))

##

fig, ax = plt.subplots(figsize=(4,3))

cmap = plu.create_palette(df_plot, 'sample', order=samples, palette='tab10')
plu.dist(df_plot, x='n_muts', by='sample', ax=ax, linewidth=.75, alpha=.3, categorical_cmap=cmap)
medians = {}
for sample in samples:
    median = np.median(df_plot.query('sample==@sample')['n_muts'])
    ax.axvline(x=median, linewidth=.75, linestyle='--', c=cmap[sample])
    medians[sample] = median

ax.text(.43, .75, str(int(medians['MDA_clones'])), transform=ax.transAxes, c=cmap['MDA_clones'])
ax.text(.38, .85, str(int(medians['MDA_lung'])), transform=ax.transAxes, c=cmap['MDA_lung'])
ax.text(.33, .95, str(int(medians['MDA_PT'])), transform=ax.transAxes, c=cmap['MDA_PT'])

plu.format_ax(ax=ax, xlabel='n enriched MT-SNVs', ylabel='Density', 
              title=f'n lentiviral clones: {df_plot.shape[0]}', reduced_spines=True)
plu.add_legend(cmap, label='Sample', loc='upper right', bbox_to_anchor=(1,1), ax=ax)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'Supp_Fig_10b.pdf'))

##