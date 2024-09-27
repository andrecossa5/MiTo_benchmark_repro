"""
Additional quality check on mito_preprocessing output. GBC library.
"""

import os
import numpy as np
import pandas as pd
from mito_utils.plotting_base import *
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro'
path_data = os.path.join(path_main, 'data', 'MI_TO_bench')
path_results = os.path.join(path_main, 'results',  'MI_TO_bench')

# Sample list
samples = ['MDA_clones', 'AML_clones', 'MDA_PT', 'MDA_lung']

# Plots
fig, axs = plt.subplots(1,4,figsize=(16,4.2))

for i,sample in enumerate(samples):
    ax = axs.ravel()[i]
    df = pd.read_csv(os.path.join(path_data, 'CBC_GBC_combinations', f'{sample}.tsv.gz'), sep='\t', index_col=0)
    sns.kdeplot(data=df, x='normalized_abundance', y='max_ratio', ax=ax)
    ax.axvline(x=.75, linestyle='--', c='r')
    ax.axhline(y=.75, linestyle='--', c='r')
    ax.set(title=sample)

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'GBC_QC', 'density_combos.png'), dpi=1000)


##


fig, axs = plt.subplots(1,4,figsize=(16,4.2))

for i,sample in enumerate(samples):
    ax = axs.ravel()[i]
    df = pd.read_csv(os.path.join(path_data, 'CBC_GBC_combinations', f'{sample}.tsv.gz'), sep='\t', index_col=0)
    scatter(df, 'umi', 'p', by='max_ratio', s=50, vmin=.2, vmax=.8, ax=ax, c='Spectral_r')
    format_ax(ax, title=sample, xlabel='nUMIs', ylabel='p')
    ax.axhline(y=1, c='k', linestyle='--')

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'GBC_QC', 'pPoisson_logumi_combos.png'), dpi=1000)


##


fig, axs = plt.subplots(1,4,figsize=(16,4.2))

for i,sample in enumerate(samples):
    ax = axs.ravel()[i]
    df = pd.read_csv(os.path.join(path_data, 'CBC_GBC_combinations', f'{sample}.tsv.gz'), sep='\t', index_col=0)
    sns.kdeplot(df.query('status=="unsupported"')['umi'], color='grey', ax=ax, fill=True, alpha=.5)
    sns.kdeplot(df.query('status=="supported"')['umi'], color='darkred', ax=ax, fill=True, alpha=.5)
    format_ax(ax, title=sample, xlabel='nUMIs', ylabel='Density')
    if i == 0:
        colors = { 'unsupported' : "grey", "supported" : 'darkred' }
        add_legend(label='CBC-GBC combo status', ax=ax, loc='upper right', bbox_to_anchor=(1,1), 
                   ticks_size=8, artists_size=8, label_size=8, colors=colors)

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'GBC_QC', 'supported_unsupported.png'), dpi=1000)



##

