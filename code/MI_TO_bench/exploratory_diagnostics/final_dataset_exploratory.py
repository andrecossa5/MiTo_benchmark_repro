"""
Additional quality check on mito_preprocessing output. MT- library.
"""

import os
import numpy as np
import pandas as pd
import random
from mito_utils.preprocessing import *
from mito_utils.plotting_base import *
from mito_utils.diagnostic_plots import *
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro'
path_data = os.path.join(path_main, 'data', 'MI_TO_bench')
path_results = os.path.join(path_main, 'results',  'MI_TO_bench')


##


# Read meta and good cells from MT-libraries
samples = ['MDA_clones', 'AML_clones', 'MDA_PT', "MDA_lung"]
meta = pd.read_csv(os.path.join(path_data, 'cells_meta.csv'), index_col=0)
L = []
for sample in samples:
    L += pd.read_csv(os.path.join(path_data, 'QC_cells', 'cells_mt', f'{sample}_mt.txt'), header=None)[0].to_list()
meta = meta.loc[meta.index.isin(L),:]

# Calculate lones frequencies
df_freq = (
    meta[['GBC', 'sample']]
    .groupby('sample')
    ['GBC'].value_counts(normalize=True)
    .reset_index(name='freq')
)


##


fig, ax = plt.subplots(figsize=(4,4))
df = meta.groupby('sample')['GBC'].nunique().to_frame('n clones').sort_values('n clones')
bar(df, 'n clones', c='grey', s=0.75, a=.7, ax=ax, edgecolor='k', annot_size=10)
format_ax(ax=ax, xticks=df.index, rotx=90, ylabel='n clones', reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'final_dataset_exploratory', 'n_clones.png'), dpi=1000)


##


fig, ax = plt.subplots(figsize=(3,4))
SH = df_freq.groupby('sample').apply(lambda x: - np.sum( x['freq'] * np.log10(x['freq'])) )
xcoor = [.9, 1.6, .4, 1]
markers = ['o', '+', 'x', 'v']
for i in range(len(samples)):
    ax.plot(xcoor[i],SH.iloc[i], marker=markers[i], c='k', alpha=.7, label=samples[i], markersize=8)
ax.set(xlim=(0,2))
ax.legend(frameon=False, bbox_to_anchor=(1,1), loc='upper left')
format_ax(ax=ax, xticks=['', 'Samples', ''], ylabel='Shannon Entropy', reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'final_dataset_exploratory', 'SH.png'), dpi=1000)


##


# Clone colors
# clones_colors = {
#     **{ clone : color for clone, color in zip(df_freq.query('sample in ["MDA_clones", "AML_clones"]')['GBC'], sc.pl.palettes.vega_20_scanpy) },
#     **{ clone : color for clone, color in zip(df_freq.query('sample in ["MDA_PT", "MDA_lung"] and freq >= 0.01')['GBC'], sc.pl.palettes.vega_20_scanpy) },
#     **{ clone : color for clone, color in zip(df_freq.query('sample in ["MDA_PT", "MDA_lung"] and freq < 0.01')['GBC'], sc.pl.palettes.zeileis_28*10) }
# }
# with open(os.path.join(path_data, 'clones_colors_sc.pickle'), 'wb') as f:
#     pickle.dump(clones_colors, f)
# 
# # Read colors
with open(os.path.join(path_data, 'clones_colors_sc.pickle'), 'rb') as f:
    clones_colors = pickle.load(f)


##


fig, axs = plt.subplots(1,4,figsize=(12,3))

order = ['AML_clones', 'MDA_clones', 'MDA_PT', 'MDA_lung']
for ax, sample in zip(axs.flat, order):
    f = .01 if sample in ['MDA_PT', 'MDA_lung'] else .001
    df_ = df_freq.query('sample==@sample and freq>=@f').set_index('GBC')
    packed_circle_plot(
        df_, covariate='freq', ax=ax, color=clones_colors, annotate=True, t_cov=.05,
        alpha=.8, linewidth=3, fontsize=7.5, fontcolor='k', fontweight='medium'
    )
    ax.set(title=sample)
    ax.axis('off')
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'circle_packed.png'), dpi=1000)


##


fig, ax = plt.subplots(figsize=(6,2.5))
meta['sample'] = pd.Categorical(meta['sample'], categories=order[::-1])
bb_plot(meta, 'sample', 'GBC', legend=False, colors=clones_colors, ax=ax)
ax.set(title='')
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'bb_plot.png'), dpi=1000)


##