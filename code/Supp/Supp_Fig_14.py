"""
Supp Fig 14
--distance_metric parameter tuning. 
"""

import os
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import mito as mt
import matplotlib.pyplot as plt
import plotting_utils as plu
import matplotlib
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_tune = os.path.join(path_main, 'data', 'bench', 'tune_distances')
path_lineage_inference = os.path.join(path_main, 'data', 'lineage_inference', 'UPMGA')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')


##


# 1. Benchmark distance metrics on a limited (n=15) number of high-quality MT-SNVs spaces ------------- # 

# Prep results list and params
RESULTS = []
max_fraction_unassigned = .1

# Iterate over samples and jobs
for sample in ['MDA_clones', 'MDA_PT', 'MDA_lung']:

    for job_id in os.listdir(os.path.join(path_lineage_inference, sample)):

        # Read AFM
        T = mt.ut.Timer()
        T.start()
        logging.info(f'Processing sample {sample}, job {job_id}...')
        afm = sc.read(os.path.join(path_lineage_inference, sample, job_id, 'afm_filtered.h5ad'))

        # Compute distances in this MT-SNV space
        del afm.obsp['distances']
        del afm.uns['distance_calculations']

        # Set metrics for comparison
        CONT_METRICS = ['euclidean', 'cosine', 'correlation']
        BIN_METRICS = ['jaccard', 'dice', 'weighted_jaccard']
        METRICS = CONT_METRICS + BIN_METRICS

        # Compute alternative cell-cell distances
        for metric in METRICS:
            mt.pp.compute_distances(afm, metric=metric, distance_key=f'{metric}_dists', precomputed=True)
        afm

        # Dists
        AUPRC = {}
        for metric in METRICS:
            labels = afm.obs['GBC'].astype('str')
            D = afm.obsp[f'{metric}_dists'].toarray()
            AUPRC[metric] = mt.ut.distance_AUPRC(D, labels)

        # kNN
        kNN_SH = {}
        kNN_PURITY = {}
        kBET = {}
        for metric in METRICS:
            labels = afm.obs['GBC'].astype('str')
            D = afm.obsp[f'{metric}_dists'].toarray()
            idx,_,_ = mt.pp.kNN_graph(D=D, k=30, from_distances=True)
            kNN_PURITY[metric] = mt.ut.NN_purity(idx, labels)
            kNN_SH[metric] = mt.ut.NN_entropy(idx, labels)
            kBET[metric] = mt.ut.kbet(idx, labels)

        # Lineage inferences and tree topology
        ARI = {}
        NMI = {}
        corr_dists = {}
        NAs = {}
        for metric in METRICS:
            labels = afm.obs['GBC'].astype('str')
            tree = mt.tl.build_tree(afm, distance_key=f'{metric}_dists', precomputed=True)
            try:
                model = mt.tl.MiToTreeAnnotator(tree)
                model.clonal_inference(max_fraction_unassigned=max_fraction_unassigned)
                test = tree.cell_meta['MiTo clone'].isna()
                inferred = tree.cell_meta.loc[~test, 'MiTo clone'].astype('str')
                GT = tree.cell_meta.loc[~test, 'GBC'].astype('str')
                NAs[metric] = test.sum() / len(test)
                ARI[metric] = mt.ut.custom_ARI(inferred, GT)
                NMI[metric] = mt.ut.normalized_mutual_info_score(inferred, GT)
            except:
                ARI[metric] = np.nan
                NMI[metric] = np.nan
                NAs[metric] = np.nan

            corr_dists[metric] = mt.ut.calculate_corr_distances(tree)[0]

        # Collect results
        df_ = pd.DataFrame({ 
            'AUPRC': AUPRC, 
            'kNN_SH': kNN_SH, 
            'kNN_PURITY': kNN_PURITY, 
            'kBET': kBET, 
            'ARI': ARI, 
            'NMI': NMI, 
            'corr_dists': corr_dists, 
            'NAs': NAs 
        })
        df_ = df_.assign(sample=sample, job_id=job_id).reset_index(names='metric')
        RESULTS.append(df_)

        logging.info(f'Finished sample {sample}, job {job_id}... in {T.stop()}.')
    

##


# Concatenate results
pd.concat(RESULTS).to_csv(os.path.join(path_tune, 'top_jobs_distances_benchmark.csv'))


##


# 2. Supp Fig 14 a,b. Evaluate different metrics ------------- # 

# Read metrics
df = pd.read_csv(os.path.join(path_tune, 'top_jobs_distances_benchmark.csv'), index_col=0)
df = df.reset_index(drop=True)

# Supp Fig 14a: Plot ARI vs distances correlation
plu.set_rcParams({'figure.dpi':100})
order = ['jaccard', 'weighted_jaccard', 'dice', 'cosine', 'correlation', 'euclidean']

fig, ax = plt.subplots(figsize=(5,3.5))

colors = plu.create_palette(df, 'metric', plu.ten_godisnot)
plu.scatter(df, 'corr_dists', 'ARI', by='metric', 
            categorical_cmap=colors, ax=ax, size=50, kwargs={'edgecolor':'k'})
plu.format_ax(ax=ax, xlabel='Distance correlation', ylabel='ARI')
plu.add_legend(colors, label='Metric', ax=ax)

fig.subplots_adjust(right=.6, bottom=.2, left=.2, top=.8)
fig.savefig(os.path.join(path_figures, 'Supp_Fig_14a.pdf'))


##


# Supp Fig 14b: all other metrics
fig, axs = plt.subplots(2,3,figsize=(5.5,3.5), sharex=True)

metrics = ['AUPRC', 'kNN_SH', 'kNN_PURITY', 'kBET', 'NMI', 'NAs']
for i, metric in enumerate(metrics):
    ax = axs.ravel()[i]
    plu.strip(df, 'metric', metric, categorical_cmap=colors, x_order=order, ax=ax)
    plu.violin(df, 'metric', metric, categorical_cmap=colors, x_order=order, ax=ax, kwargs={'alpha':.5})
    plu.format_ax(ax=ax, rotx=90, xlabel='', reduced_spines=True)

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'Supp_Fig_14b.pdf'))


##


# 3. Supp Fig 14c. Visualize example cases ------------- # 

samples = ['MDA_PT', 'MDA_clones']
job_ids = ['39aaaaf933','b9f2160731']

for sample, job_id in zip(samples, job_ids):

    # Read
    afm = sc.read(os.path.join(path_lineage_inference, sample, job_id, 'afm_filtered.h5ad'))

    # Compute distances in this MT-SNV space
    del afm.obsp['distances']
    del afm.uns['distance_calculations']

    # Compute alternative cell-cell distances
    order = ['jaccard', 'weighted_jaccard', 'dice', 'cosine', 'correlation', 'euclidean']
    for metric in order:
        mt.pp.compute_distances(afm, metric=metric, distance_key=f'{metric}_dists', precomputed=True)

    # Cell phylogeny
    fig, axs = plt.subplots(1,6,figsize=(12,2.2))

    for i,x in enumerate(order):
        ax = axs.ravel()[i]
        tree = mt.tl.build_tree(afm, precomputed=True, distance_key=f'{x}_dists')
        mt.pl.heatmap_distances(afm, distance_key=f'{x}_dists', tree=tree, ax=ax)
        ax.set(title=x)

    fig.tight_layout()
    fig.savefig(os.path.join(path_figures, f'Supp_Fig_14_c_{sample}.pdf'))


##


# 3. Supp Fig 14c. Plot metrics with standard or weighted jaccard distance --------------------- # 

L = []
for folder,_,files in os.walk(path_tune):
    if any([ x.startswith('all') for x in files]):
        df,metrics,options = mt.ut.format_tuning(folder)
        L.append(df)

df = pd.concat(L)

varying_options = (df[options].nunique()).loc[lambda x:x>1].index.to_list()
metrics_of_interest = ['ARI', 'NMI', 'corr', 'AUPRC', 'mean_CI', 'mean_RI', 'n MiTo clone']
# (
#     df.groupby(['sample']+varying_options)
#     [['ARI', 'NMI', 'corr', 'AUPRC', 'n_cells', 'n_vars', 'n_GBC_groups']]
#     .median()
# )

##


# Plot
fig, axs =  plt.subplots(1,len(metrics_of_interest),figsize=(14,5))

x_order = ['MDA_clones', 'MDA_lung', 'MDA_PT']
by_order = ['jaccard', 'weighted_jaccard']

cmap = plu.create_palette(df, 'metric', order=by_order, palette='Reds')
for i,metric in enumerate(metrics_of_interest):
    plu.bar(
        df.query('min_AD=="2"'), 
        x='sample', y=metric, 
        by='metric',
        x_order=x_order,
        by_order=by_order,
        categorical_cmap=cmap,
        ax=axs[i]
    )
    plu.format_ax(ax=axs[i], xlabel='', ylabel=metric, reduced_spines=True, rotx=90)

plu.add_legend(cmap, ax=axs[3], label='Distance metric', loc='center', bbox_to_anchor=(.5,1.2), ncols=2)
fig.subplots_adjust(top=.75, left=.1, bottom=.25, right=.9, wspace=.6)
fig.savefig(os.path.join(path_figures, 'Supp_Fig_14d.pdf'))


##