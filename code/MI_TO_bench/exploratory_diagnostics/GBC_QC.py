"""
Additional quality check on mito_preprocessing output. GBC library.
"""

import os
import pandas as pd
from mito_utils.plotting_base import *


##


# Set plotting options
set_rcParams()
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro'
path_data = os.path.join(path_main, 'data', 'MI_TO_bench')
path_results = os.path.join(path_main, 'results',  'MI_TO_bench')

# Sample list
samples = ['MDA_clones', 'MDA_PT', 'MDA_lung']

# Plots
fig, axs = plt.subplots(2,3,figsize=(9,6))

for i,sample in enumerate(samples):
    ax = axs[0,i]
    df = pd.read_csv(os.path.join(path_data, 'CBC_GBC_combinations', f'{sample}.tsv.gz'), sep='\t', index_col=0)
    sns.kdeplot(data=df, x='normalized_abundance', y='max_ratio', ax=ax)
    ax.axvline(x=.75, linestyle='--', c='r')
    ax.axhline(y=.75, linestyle='--', c='r')
    ax.set_title(sample)                           # Use default font size from rcParams
    ax.set_xlabel('Normalized abundance')
    ax.set_ylabel('Max ratio')

for i,sample in enumerate(samples):
    ax = axs[1,i]
    df = pd.read_csv(os.path.join(path_data, 'CBC_GBC_combinations', f'{sample}.tsv.gz'), sep='\t', index_col=0)
    sns.kdeplot(df.query('status=="unsupported"')['umi'], color='grey', ax=ax, fill=True, alpha=.5)
    sns.kdeplot(df.query('status=="supported"')['umi'], color='darkred', ax=ax, fill=True, alpha=.5)
    format_ax(ax=ax, xlabel='nUMIs', ylabel='Density', xlabel_size=10, ylabel_size=10)

    # Add legend only to the first plot in the second row
    if i == 0:
        colors = { 'unsupported' : "grey", "supported" : 'darkred' }
        add_legend(label='CB-GBC combo status', ax=ax, loc='upper right', bbox_to_anchor=(1,1), 
                   ticks_size=10, artists_size=10, label_size=10, colors=colors)


fig.tight_layout()
fig.savefig(os.path.join(path_results, 'GBC_QC', 'GBC_QC.pdf'))


##


fig, axs = plt.subplots(1,3,figsize=(12,4.2))

for i,sample in enumerate(samples):
    ax = axs.ravel()[i]
    df = pd.read_csv(os.path.join(path_data, 'CBC_GBC_combinations', f'{sample}.tsv.gz'), sep='\t', index_col=0)
    scatter(df, 'umi', 'p', by='max_ratio', s=50, vmin=.2, vmax=.8, ax=ax, c='Spectral_r')
    format_ax(ax, title=sample, xlabel='nUMIs', ylabel='p', xlabel_size=15, title_size=15, ylabel_size=15)
    ax.axhline(y=1, c='k', linestyle='--')

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'GBC_QC', 'pPoisson_logumi_combos.pdf'))



##