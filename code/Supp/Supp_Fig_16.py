
"""
Supp Fig 16
Clonal inference methods benchmaring (nf-MiTo BENCH).
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import mito as mt
import matplotlib
import matplotlib.pyplot as plt
import plotting_utils as plu
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_bench = os.path.join(path_main, 'data', 'bench', 'clonal_inference')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')

# Set visualization params
plu.set_rcParams({'figure.dpi':350})


##


# Utils
def _format_time(s):

    time = None
    if any([x=='m' for x in s]):
        minutes, seconds = s.split(' ')
        minutes = float(minutes[:-1]) * 60
        seconds = float(seconds[:-1])
        time = minutes + seconds
    else:
        time = float(s[:-1])

    return time

##

def _format_memory(s):

    try:
        memory, unit = s.split(' ')
        if unit=='MB':
            memory = float(memory) * .001
        else:
            memory = float(memory)
    except:
        memory = np.nan

    return memory

##


# Benchmark time and memory usage across clonal inference methods --------------------------------------- # 

# Read traces
samples = ['MDA_clones', 'MDA_lung', 'MDA_PT']

L = []
for sample in samples:
    df = pd.read_csv(os.path.join(path_bench, f'{sample}.txt'), sep='\t', index_col=0)  
    L.append(df)

df = pd.concat(L)
df = df[['name', 'duration', 'realtime', 'peak_rss', '%cpu']].copy()

# Format values
mapping = {'LEIDEN':'leiden', 'VIREO':'vireoSNP', 'MITO_BENCH':'MiTo', 'CCLONE':'CClone'}
df['method'] = df['name'].map(lambda x: x.split(' ')[0].split(':')[-1]).map(mapping)
df['duration'] = df['duration'].map(lambda x: _format_time(x)) / 60
df['peak_rss'] = df['peak_rss'].map(lambda x: _format_memory(x)) / 1000
df['sample'] = df['name'].map(lambda x: x.split(' ')[1].replace('(','').replace(':', ''))
df = df.set_index('name')


##


# Viz
df_bench = mt.ut.extract_bench_df(path_bench)
order = df_bench.groupby('method')['ARI'].median().sort_values().index
cmap = plu.create_palette(df, 'method', order=order, col_list=sc.pl.palettes.vega_10_scanpy)

# Fig
fig, axs = plt.subplots(1,2, figsize=(6,3))

plu.strip(df, x='sample', y='duration', by='method', x_order=samples, categorical_cmap=cmap, by_order=order, ax=axs[0])
plu.format_ax(ax=axs[0], ylabel='Time (min)', xlabel='', reduced_spines=True)
plu.add_legend(cmap, label='Method', ax=axs[0], bbox_to_anchor=(0,1))
plu.strip(df, x='sample', y='peak_rss', by='method', x_order=samples, categorical_cmap=cmap, by_order=order, ax=axs[1])
plu.format_ax(ax=axs[1], ylabel='Peak rss (GB)', xlabel='', reduced_spines=True)

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'Supp_Fig_16.pdf'))


##