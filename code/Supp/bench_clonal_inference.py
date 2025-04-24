"""
Bench clonal inference.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import pickle
import mito as mt
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import plotting_utils as plu
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_figures = os.path.join(path_main, 'results', 'figures', 'Fig3')
path_results = os.path.join(path_main, 'results', 'others', 'Fig3')


##


# 1. Visualize MT-SNVs from top chosen MT-SNVs spaces --------------------------------------- # 

# Params
plu.set_rcParams()
matplotlib.rcParams.update({'figure.dpi':150})
path_data = os.path.join(path_main, 'data', 'lineage_inference', 'UPMGA')

# Samples order
samples = ['MDA_clones', 'MDA_PT', 'MDA_lung']

# Here we go
fig, axs = plt.subplots(1,3,figsize=(5,2))

for i,sample in enumerate(samples):
    
    n_positive = []
    mean_AD_in_positives = []
    variants = []

    for job in os.listdir(os.path.join(path_data, sample)):
        afm = sc.read(os.path.join(path_data, sample, job, f'afm_filtered.h5ad'))
        variants += afm.var_names.to_list()
        n_positive += (afm.layers['bin'].toarray()==1).sum(axis=0).tolist()
        mean_AD_in_positives += np.nanmean(
            np.where(afm.layers['bin'].toarray()==1, afm.layers['AD'].toarray(), np.nan), axis=0
        ).tolist()

    df = (
        pd.DataFrame(
            {'n_positive':n_positive, 'mean_AD_in_positives':mean_AD_in_positives}, 
            index=variants
        )
        .reset_index(names='var')
    )

    ax = axs.ravel()[i]
    sns.kdeplot(data=df, x='n_positive', y='mean_AD_in_positives', fill=False, color='#41767F', ax=ax)
    sns.kdeplot(data=df, x='n_positive', y='mean_AD_in_positives', fill=True, color='#41767F', alpha=.7, ax=ax)
    n_vars = df['var'].unique().size
    median_n_positives = df['n_positive'].median()
    median_mean_AD_in_positives = df['mean_AD_in_positives'].median()
    plu.format_ax(
        title=f'{sample}: {n_vars}', 
        xlabel='n +cells', 
        ylabel='Mean ALT\nUMIs +cells' if i==0 else '', 
        ax=ax, 
        reduced_spines=True
    )
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    ax.vlines(x=median_n_positives, ymin=y_min, ymax=median_mean_AD_in_positives, colors='red', linestyles='dashed')
    ax.hlines(y=median_mean_AD_in_positives, xmin=x_min, xmax=median_n_positives, colors='red', linestyles='dashed')

    ax.plot(median_n_positives, median_mean_AD_in_positives, 'rx', markersize=10)

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'selected_MT-SNVs.pdf'))


##


# 2. Visualize ARI, NMI across samples and clonal inference methods --------------------------------------- # 

# Extract bench dataframe
path_data = os.path.join(path_main, 'data', 'bench', 'clonal_inference')
df_bench = mt.ut.extract_bench_df(path_data)
df_bench.groupby(['sample', 'method'])[['ARI', 'NMI']].median()


##


# Viz performance

# Fig
fig, axs = plt.subplots(1,2,figsize=(5.5,2.5))

df_bench['sample'] = pd.Categorical(
    df_bench['sample'], categories=['MDA_clones', 'MDA_lung', 'MDA_PT']
)
order = df_bench.groupby('method')['ARI'].median().sort_values().index
colors = plu.create_palette(df, 'method', order=order, col_list=sc.pl.palettes.vega_10_scanpy)

# ARI
plu.strip(df_bench, 'sample', 'ARI', by='method', by_order=order, categorical_cmap=colors, ax=axs[0])
plu.bar(df_bench, 'sample', 'ARI', by='method', by_order=order, categorical_cmap=colors, ax=axs[0])
axs[0].set_ylim((-.02,1))
axs[0].axhline(.9, linestyle='--', color='k', linewidth=.5)
plu.format_ax(ax=axs[0], xlabel='', ylabel='ARI', reduced_spines=True, rotx=90)

# NMI
plu.strip(df_bench, 'sample', 'NMI', by='method', by_order=order, categorical_cmap=colors, ax=axs[1])
plu.bar(df_bench, 'sample', 'NMI', by='method', by_order=order, categorical_cmap=colors, ax=axs[1])
axs[1].set_ylim((-.02,1))
axs[1].axhline(.9, linestyle='--', color='k', linewidth=.5)
plu.format_ax(ax=axs[1], xlabel='', ylabel='NMI', reduced_spines=True, rotx=90)
plu.add_legend(colors, label='Method', ax=axs[1])

# Readjust and save
fig.subplots_adjust(right=.75, top=.9, left=.1, bottom=.35, wspace=.4)
fig.savefig(os.path.join(path_figures, 'clonal_reconstruction_performance.pdf'))


##


# 3. Visualize clonal assignment across samples and clonal inference methods --------------------------------------- # 

# Set paths
path_lineage_inference = os.path.join(path_main, 'data', 'lineage_inference', 'UPMGA')
path_bench = os.path.join(path_main, 'data', 'bench', 'clonal_inference')

# Extract bench dataframe
df_bench = mt.ut.extract_bench_df(path_bench)
df_bench.groupby(['sample', 'method']).describe()

# Choose one job per sample, for visualization purposes
top_3_jobs = (
    df_bench
    .groupby(['sample', 'job_id'])
    ['ARI'].median().sort_values().reset_index()
    .groupby('sample')
    .apply(lambda x: x['job_id'].values[x['ARI'].argmax()])
    .to_dict()
)
top_3_jobs

# Load colors
path_colors = os.path.join(path_main, 'data', 'general')
with open(os.path.join(path_colors, 'clones_colors_sc.pickle'), 'rb') as f:
    clone_colors = pickle.load(f)


##


# Viz
for i,sample in enumerate(['MDA_clones', 'MDA_PT', 'MDA_lung']):

    # job_id
    job_id = top_3_jobs[sample]

    # Read
    afm = sc.read(os.path.join(path_lineage_inference, sample, job_id, 'afm_filtered.h5ad'))
    mt.pp.reduce_dimensions(afm, metric='cosine')

    n_cells = afm.shape[0]
    n_vars = afm.shape[1]
    n_GBC = afm.obs['GBC'].unique().size
    print(f'{sample}: {n_cells}, {n_vars}, {n_GBC}')

    # Ground Truth: UMAP
    f,ax = plt.subplots(figsize=(4,4))
    s = None if sample != 'MDA_clones' else 200
    mt.pl.draw_embedding(afm, feature='GBC', ax=ax, categorical_cmap=clone_colors, size=s)
    f.tight_layout()
    f.savefig(os.path.join(path_figures, f'{sample}_GT_UMAP.pdf'))
    
    # Confusion matrices
    fig, axs = plt.subplots(1,4,figsize=(8.5,2.2))

    for j,method in enumerate(['MiTo', 'vireoSNP', 'leiden', 'CClone']):

        with open(os.path.join(path_bench, sample, f'{job_id}_{method}.pickle'), 'rb') as f:
            d = pickle.load(f)

        labels = d['labels']
        afm.obs[method] = labels
        afm.obs.loc[afm.obs[method].isna(), method] = 'unassigned'
        ARI = d['ARI']
        NMI = d['NMI']

        df_plot = pd.crosstab(afm.obs['GBC'], afm.obs[method])
        order_row = afm.obs['GBC'].value_counts().index
        order_col = afm.obs.loc[afm.obs[method]!='unassigned', method].value_counts().index
        df_plot = df_plot.loc[order_row,order_col]
        plu.plot_heatmap(
            df_plot, ax=axs[j], palette='Blues', label='n cells', 
            x_names=False, y_names=False, xlabel='Inferred', ylabel='Ground Truth',
            title=f'{method} (n clones={df_plot.shape[1]})\nARI: {ARI:.2f}, NMI: {NMI:.2f}'
        )
 
        fig.tight_layout()
        fig.savefig(os.path.join(path_figures, f'{sample}_confusion.pdf'))


##


# 4. Benchmark time and memory --------------------------------------- # 

samples = ['MDA_clones', 'MDA_lung', 'MDA_PT']
path_bench = os.path.join(path_main, 'data', 'bench', 'clonal_inference')

# Read Nextflow traces
L = []
for sample in samples:
    df = pd.read_csv(os.path.join(path_bench, f'{sample}.txt'), sep='\t', index_col=0)  
    L.append(df)

df = pd.concat(L)
df = df[['name', 'duration', 'realtime', 'peak_rss', '%cpu']].copy()


##

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


# Format values
mapping = {'LEIDEN':'leiden', 'VIREO':'vireoSNP', 'MITO_BENCH':'MiTo', 'CCLONE':'CClone'}
df['method'] = df['name'].map(lambda x: x.split(' ')[0].split(':')[-1]).map(mapping)
df['duration'] = df['duration'].map(lambda x: _format_time(x)) / 60
df['peak_rss'] = df['peak_rss'].map(lambda x: _format_memory(x)) / 1000
df['sample'] = df['name'].map(lambda x: x.split(' ')[1].replace('(','').replace(':', ''))
df = df.set_index('name')


##


# Viz

# Colors
df_bench = mt.ut.extract_bench_df(path_bench)
order = df_bench.groupby('method')['ARI'].median().sort_values().index
cmap = plu.create_palette(df, 'method', order=order, col_list=sc.pl.palettes.vega_10_scanpy)

fig, axs = plt.subplots(1,2, figsize=(6,3))

plu.strip(df, x='sample', y='duration', by='method', x_order=samples,
          categorical_cmap=cmap, by_order=order, ax=axs[0])
plu.format_ax(ax=axs[0], ylabel='Time (min)', xlabel='', reduced_spines=True)
plu.add_legend(cmap, label='Method', ax=axs[0], bbox_to_anchor=(0,1))

plu.strip(df, x='sample', y='peak_rss', by='method', x_order=samples,
          categorical_cmap=cmap, by_order=order, ax=axs[1])
plu.format_ax(ax=axs[1], ylabel='Peak rss (GB)', xlabel='', reduced_spines=True)


fig.tight_layout()
plt.show()


##