"""
Fig.3
Clonal inference benchmark.
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
# matplotlib.use('macOSX')          # On macOS only


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_bench = os.path.join(path_main, 'data', 'bench', 'clonal_inference') 
path_filtered_afms = os.path.join(path_main, 'data', 'lineage_inference', 'UPMGA') 
path_figures = os.path.join(path_main, 'results', 'figures', 'Fig3')
path_results = os.path.join(path_main, 'results', 'others', 'Fig3')


##


# Fig 3a. Visualize MT-SNVs from MT-SNVs spaces chosen for benchmarking ---------------------------- # 

# Params
plu.set_rcParams({'figure.dpi':350})

# Samples order
samples = ['MDA_clones', 'MDA_PT', 'MDA_lung']

# Plot
fig, axs = plt.subplots(1,3,figsize=(8,3))

for i,sample in enumerate(samples):
    
    n_positive = []
    mean_AD_in_positives = []
    variants = []

    for job in os.listdir(os.path.join(path_filtered_afms, sample)):
        afm = sc.read(os.path.join(path_filtered_afms, sample, job, 'afm_filtered.h5ad'))
        variants += afm.var_names.to_list()
        n_positive += (afm.layers['bin']==1).sum(axis=0).A1.tolist()
        mean_AD_in_positives += (
            np.nanmean(
                np.where(afm.layers['bin'].toarray()==1, afm.layers['AD'].toarray(), np.nan), 
                axis=0
            )
            .tolist()
        )

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
fig.savefig(os.path.join(path_figures, 'Fig_3a.pdf'))


##


# Fig 2b. Visualize key metrics samples and clonal inference methods --------------------------------------- # 

# Extract bench dataframe
df_bench = mt.ut.extract_bench_df(path_bench)
df_bench.groupby('method').describe().T

# Add 'nMT - nGBC' column to df_bench
results = []
for i,sample in enumerate(samples):
    for job in os.listdir(os.path.join(path_filtered_afms, sample)):
        afm = sc.read(os.path.join(path_filtered_afms, sample, job, 'afm_filtered.h5ad'))
        results.append([job, afm.obs['GBC'].nunique()])

df_bench = df_bench.merge(
    pd.DataFrame(results, columns=['job_id', 'n_GBC']), 
    on='job_id', how='left'
)
df_bench['nMT - nGBC'] = df_bench['n_inferred'] - df_bench['n_GBC']


##


# Plot
fig, axs = plt.subplots(1,2,figsize=(9,4))

df_bench['sample'] = pd.Categorical(
    df_bench['sample'], categories=['MDA_clones', 'MDA_lung', 'MDA_PT']
)
order = df_bench.groupby('method')['ARI'].median().sort_values().index
colors = { k:v for k,v in zip(order, sc.pl.palettes.vega_10_scanpy) }

# 'nMT - nGBC'
plu.strip(df_bench, 'sample', 'nMT - nGBC', by='method', by_order=order, categorical_cmap=colors, ax=axs[0])
plu.bar(df_bench, 'sample', 'nMT - nGBC', by='method', by_order=order, categorical_cmap=colors, ax=axs[0])
plu.format_ax(ax=axs[0], xlabel='', ylabel='nMT - nGBC', reduced_spines=True)

# ARI
plu.strip(df_bench, 'sample', 'ARI', by='method', by_order=order, categorical_cmap=colors, ax=axs[1])
plu.bar(df_bench, 'sample', 'ARI', by='method', by_order=order, categorical_cmap=colors, ax=axs[1])
axs[1].set_ylim((-.02,1))
axs[1].axhline(.9, linestyle='--', color='k', linewidth=.5)
plu.format_ax(ax=axs[1], xlabel='', ylabel='ARI', reduced_spines=True)

fig.subplots_adjust(right=.75, top=.85, left=.15, bottom=.15)
fig.savefig(os.path.join(path_figures, 'Fig_3b.pdf'))


##


# Fig 3c. Visualize clonal assignment across samples and clonal inference methods --------------------------------------- # 

# Extract bench dataframe
df_bench = mt.ut.extract_bench_df(path_bench)

# Choose one representative job per sample (visualization purposes)
top_3_jobs = (
    df_bench
    .groupby(['sample', 'job_id'])
    ['ARI'].median().sort_values().reset_index()
    .groupby('sample')
    .apply(lambda x: x['job_id'].values[x['ARI'].argmax()])
    .to_dict()
)

# Load clones colors
path_colors = os.path.join(path_main, 'data', 'general')
with open(os.path.join(path_colors, 'clones_colors_sc.pickle'), 'rb') as f:
    clone_colors = pickle.load(f)


##


# Plot
samples = ['MDA_clones', 'MDA_PT', 'MDA_lung']
methods = ['MiTo', 'vireoSNP', 'leiden', 'CClone']

for i,sample in enumerate(samples):

    # job_id
    job_id = top_3_jobs[sample]

    # Read AFM
    afm = sc.read(os.path.join(path_filtered_afms, sample, job_id, 'afm_filtered.h5ad'))
    mt.pp.reduce_dimensions(afm)

    # UMAP
    fig, ax = plt.subplots(figsize=(4,4))
    s = None if sample != 'MDA_clones' else 200
    mt.pl.draw_embedding(afm, feature='GBC', ax=ax, categorical_cmap=clone_colors, size=s)
    fig.tight_layout()
    fig.savefig(os.path.join(path_figures, f'Fig_3c_{sample}_UMAP.pdf'))
    
    # Confusion matrices
    fig, axs = plt.subplots(1,4,figsize=(12,3))
    
    # Add MT-clones labels from each methods to afm.obs
    for j,method in enumerate(methods):

        # Load labels
        with open(os.path.join(path_bench, sample, f'{job_id}_{method}.pickle'), 'rb') as f:
            d = pickle.load(f)
        afm.obs[method] = d['labels']
        afm = afm[~afm.obs[method].isna()].copy()

        # Confusion matrix
        df_plot = pd.crosstab(afm.obs['GBC'].astype('str'), afm.obs[method].astype('str'))
        order_row = afm.obs['GBC'].value_counts().index
        order_col = afm.obs[method].value_counts().index
        df_plot = df_plot.loc[order_row,order_col]

        # Plot
        ARI, NMI = (
            df_bench.query('sample==@sample and method==@method and job_id==@job_id')
            [['ARI', 'NMI']].values[0]
        )
        plu.plot_heatmap(
            df_plot, ax=axs[j], palette='Blues', label='n cells',
            x_names=False, y_names=False, xlabel='Inferred', ylabel='Ground Truth',
            title=f'{method} (n clones={afm.obs[method].nunique()})\nARI: {ARI:.2f}, NMI: {NMI:.2f}'
        )

    fig.tight_layout()
    fig.savefig(os.path.join(path_figures, f'Fig_3c_{sample}_confusion.pdf'))


##