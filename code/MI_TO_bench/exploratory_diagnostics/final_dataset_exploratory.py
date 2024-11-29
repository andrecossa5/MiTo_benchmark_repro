"""
Visualizatio of the final, QCed multi-modal MiTo benchmarking dataset.
"""

import os
import numpy as np
import pandas as pd
from mito_utils.preprocessing import *
from mito_utils.plotting_base import *
from mito_utils.diagnostic_plots import *
matplotlib.use('macOSX')
set_rcParams()


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro'
path_data = os.path.join(path_main, 'data', 'MI_TO_bench')
path_results = os.path.join(path_main, 'results',  'MI_TO_bench')


##


# Samples
samples = ['MDA_clones', 'MDA_PT', 'MDA_lung']

# Read meta
meta = pd.read_csv(os.path.join(path_data, 'cells_meta.csv'), index_col=0)
L = []
for sample in samples:
    L += pd.read_csv(os.path.join(path_data, 'QC_cells', f'{sample}.txt'), header=None)[0].to_list()
meta = meta.loc[meta.index.isin(L),:]
meta['sample'] = meta['sample'].astype('str')
meta['GBC'] = meta['GBC'].astype('str')

# Calculate clones frequencies
df_freq = (
    meta[['GBC', 'sample']]
    .groupby('sample')
    ['GBC'].value_counts(normalize=True)
    .reset_index(name='freq')
)


##


fig, axs = plt.subplots(1,4,figsize=(9,3.5))

fontsize = {'xlabel_size':10, 'ylabel_size':10, 'xticks_size':10}

ax = axs[0]
df = meta.groupby('sample').size().to_frame('n cells').sort_values('n cells')
bar(df, 'n cells', c='#B7B4B0', s=0.75, a=.7, ax=ax, edgecolor='k', annot_size=10)
format_ax(ax=ax, xticks=df.index, rotx=90, ylabel='n cells', reduced_spines=True, **fontsize)

ax = axs[1]
df = meta.groupby('sample')['GBC'].nunique().to_frame('n clones').sort_values('n clones')
bar(df, 'n clones', c='#B7B4B0', s=0.75, a=.7, ax=ax, edgecolor='k', annot_size=10)
format_ax(ax=ax, xticks=df.index, rotx=90, ylabel='n clones', reduced_spines=True)

ax = axs[2]
df = meta.groupby(['sample','GBC']).size().to_frame('n cells').reset_index().sort_values('n cells')
df['n cells'] = np.log10(df['n cells'])
box(df, 'sample', 'n cells', c='white', ax=ax, order=['MDA_clones', 'MDA_lung', 'MDA_PT'])
format_ax(ax=ax, rotx=90, ylabel='log(n cells)', reduced_spines=True, **fontsize)

ax = axs[3]
SH = df_freq.groupby('sample').apply(lambda x: - np.sum( x['freq'] * np.log2(x['freq'])) ).loc[samples]
xcoor = [.7, 1, 1]
markers = ['o', '+', 'x']
for i in range(len(samples)):
    ax.plot(xcoor[i],SH.iloc[i], marker=markers[i], c='k', alpha=.7, label=samples[i], markersize=8)
ax.set(xlim=(0,2))
ax.legend(frameon=False, bbox_to_anchor=(1,1), loc='upper left')
format_ax(ax=ax, xticks=['', 'Samples', ''], ylabel='Shannon Entropy', reduced_spines=True)

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'final_dataset_exploratory', 'cells_and_clones.pdf'))


##


# Clone colors
# clones_colors = {
#     **{ clone : color for clone, color in zip(df_freq.query('sample in ["MDA_clones", "AML_clones"]')['GBC'], sc.pl.palettes.vega_20_scanpy) },
#     **{ clone : color for clone, color in zip(df_freq.query('sample in ["MDA_PT", "MDA_lung"] and freq >= 0.01')['GBC'], sc.pl.palettes.vega_20_scanpy) },
#     **{ clone : color for clone, color in zip(df_freq.query('sample in ["MDA_PT", "MDA_lung"] and freq < 0.01')['GBC'], sc.pl.palettes.zeileis_28*10) }
# }
# with open(os.path.join(path_data, 'clones_colors_sc.pickle'), 'wb') as f:
#     pickle.dump(clones_colors, f)

# Read colors
with open(os.path.join(path_data, 'clones_colors_sc.pickle'), 'rb') as f:
    clones_colors = pickle.load(f)


##


fig, axs = plt.subplots(1,3,figsize=(10,3))

order = ['MDA_clones', 'MDA_PT', 'MDA_lung']
for ax, sample in zip(axs.flat, order):
    f = .01 if sample in ['MDA_PT', 'MDA_lung'] else .001
    df_ = df_freq.query('sample==@sample and freq>=@f').set_index('GBC')
    packed_circle_plot(
        df_, covariate='freq', ax=ax, color=clones_colors, annotate=True, t_cov=.05,
        alpha=.8, linewidth=2, fontsize=8.5, fontcolor='k', fontweight='medium'
    )
    ax.set(title=sample)
    ax.axis('off')

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'final_dataset_exploratory', 'circle_packed.pdf'))


##


fig, ax = plt.subplots(figsize=(6,2.5))
meta['sample'] = pd.Categorical(meta['sample'], categories=order[::-1])
bb_plot(meta, 'sample', 'GBC', legend=False, colors=clones_colors, ax=ax)
ax.set(title='')
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'final_dataset_exploratory', 'bb_plot.pdf'))


##