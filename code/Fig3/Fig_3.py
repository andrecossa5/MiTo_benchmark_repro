"""
Main Fig.3. Benchmarking clonal inference. MiTo vs clonal inference methods

1. Chosen MT-SNVs range. 
2. ARI,NMI
3. Clonal assignments (UMAPs/confusion matrices)
"""

import os
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import mito as mt
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import plotting_utils as plu
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'data', 'bench', 'clonal_inference') 
path_figures = os.path.join(path_main, 'results', 'figures', 'Fig3')
path_results = os.path.join(path_main, 'results', 'others', 'Fig3')


##


# Utils  --------------------------------------- # 

def extract_bench_df(path_data):
    """
    Extract performances from benchmarking folder.
    """
    bench_results = []
    for folder, _, files in os.walk(path_data):
        for file in files:
            splitted = os.path.join(folder, file).split('/')
            sample = splitted[-3]
            job_id = splitted[-2]
            if file.startswith('bench'):
                with open(os.path.join(folder, file), 'rb') as f:
                    d = pickle.load(f)
                res = {}
                res['ARI'] = [ d[k]['ARI'] for k in d ]
                res['NMI'] = [ d[k]['NMI'] for k in d ]
                res['% unassigned'] = [ d[k]['% unassigned'] for k in d ]
                res['sample'] = [ sample for _ in range(len(d)) ]
                res['job_id'] = [ job_id for _ in range(len(d)) ]
                res['method'] = list(d.keys())
                bench_results.append(pd.DataFrame(res))

    df_bench = pd.concat(bench_results).reset_index()

    return df_bench


##


# 1. Visualize MT-SNVs from top chosen MT-SNVs spaces --------------------------------------- # 

# Params
plu.set_rcParams()

# Samples order
samples = ['MDA_clones', 'MDA_PT', 'MDA_lung']

# Here we go
fig, axs = plt.subplots(1,3,figsize=(10,4))

for i,sample in enumerate(samples):
    
    n_positive = []
    mean_AD_in_positives = []
    variants = []

    for job in os.listdir(os.path.join(path_data, sample)):
        afm = sc.read(os.path.join(path_data, sample, job, 'afm.h5ad'))
        variants += afm.var_names.to_list()
        n_positive += (afm.layers['bin'].A==1).sum(axis=0).tolist()
        mean_AD_in_positives += np.nanmean(
            np.where(afm.layers['bin'].A==1, afm.layers['AD'].A, np.nan), axis=0
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
    plu.format_ax(title=f'{sample}: {n_vars} MT-SNVs', 
                  xlabel='n +cells', ylabel='Mean n ALT UMIs in +cells', ax=ax, reduced_spines=True)
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
df_bench = extract_bench_df(path_data)
df_bench.groupby(['sample', 'method']).describe()


##


# Viz performance

# Fig
fig, axs = plt.subplots(1,2,figsize=(9,4))

df_bench['sample'] = pd.Categorical(
    df_bench['sample'], categories=['MDA_clones', 'MDA_lung', 'MDA_PT']
)
order = df_bench.groupby('method')['ARI'].median().sort_values().index
colors = { k:v for k,v in zip(order, sc.pl.palettes.vega_10_scanpy) }

# ARI
plu.strip(df_bench, 'sample', 'ARI', by='method', by_order=order, categorical_cmap=colors, ax=axs[0])
plu.bar(df_bench, 'sample', 'ARI', by='method', by_order=order, categorical_cmap=colors, ax=axs[0])
axs[0].set_ylim((-.02,1))
axs[0].axhline(.9, linestyle='--', color='k', linewidth=.5)
plu.format_ax(ax=axs[0], xlabel='', ylabel='ARI', reduced_spines=True)

# NMI
plu.strip(df_bench, 'sample', 'NMI', by='method', by_order=order, categorical_cmap=colors, ax=axs[1])
plu.bar(df_bench, 'sample', 'NMI', by='method', by_order=order, categorical_cmap=colors, ax=axs[1])
axs[1].set_ylim((-.02,1))
axs[1].axhline(.9, linestyle='--', color='k', linewidth=.5)
plu.format_ax(ax=axs[1], xlabel='', ylabel='NMI', reduced_spines=True)

# Readjust and save
fig.subplots_adjust(right=.75, top=.85, left=.15, bottom=.15)
fig.savefig(os.path.join(path_figures, 'clonal_reconstruction_performance.pdf'))


##


# 3. Visualize clonal assignment across samples and clonal inference methods --------------------------------------- # 

# Extract bench dataframe
df_bench = extract_bench_df(path_data)
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
    afm = sc.read(os.path.join(path_data, sample, job_id, 'afm.h5ad'))
    mt.pp.reduce_dimensions(afm)

    # Ground Truth: UMAP
    f,ax = plt.subplots(figsize=(4,4))
    s = None if sample != 'MDA_clones' else 200
    mt.pl.draw_embedding(afm, feature='GBC', ax=ax, categorical_cmap=clone_colors, size=s)
    f.tight_layout()
    f.savefig(os.path.join(path_figures, f'{sample}_GT_UMAP.pdf'))
    
    # Confusion matrices
    fig, axs = plt.subplots(1,4,figsize=(12,3))

    with open(os.path.join(path_data, sample, job_id, 'bench_clonal_recontruction.pickle'), 'rb') as f:
        d = pickle.load(f)
    
    for j,method in enumerate(['MiTo', 'vireoSNP', 'leiden', 'CClone']):

        labels = d[method]['labels']
        afm.obs[method] = labels
        afm.obs[method][afm.obs[method].isna()] = 'unassigned'
        afm.obs[method] = pd.Categorical(labels)
        ARI, NMI = (
            df_bench.query('sample==@sample and method==@method and job_id==@job_id')
            [['ARI', 'NMI']].values[0]
        )

        df_plot = pd.crosstab(afm.obs['GBC'].astype('str'), afm.obs[method].astype('str'))
        order_row = afm.obs.loc[afm.obs['GBC'].isin(df_plot.index)]['GBC'].value_counts().index
        order_col = afm.obs.loc[afm.obs[method].isin(df_plot.columns)][method].value_counts().index
        df_plot = df_plot.loc[order_row,order_col]
        plu.plot_heatmap(
            df_plot, ax=axs[j], palette='Blues', label='n cells',
            x_names=False, y_names=False, xlabel='Inferred', ylabel='Ground Truth',
            title=f'Inferred clones: {df_plot.shape[1]}\nARI: {ARI:.2f}, NMI: {NMI:.2f}'
        )

        fig.tight_layout()
        fig.savefig(os.path.join(path_figures, f'{sample}_confusion.pdf'))


##