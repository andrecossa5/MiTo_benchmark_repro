"""
Bench pp_method
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
path_data = os.path.join(path_main, 'data', 'bench', 'tune_pp_method')
# path_figures = os.path.join(path_main, 'results', 'figures', 'Fig2')
# path_results = os.path.join(path_main, 'results', 'others', 'Fig2')


##


# Format metrics df 
L = []
for folder,_,files in os.walk(path_data):
    if any([ x.startswith('all') for x in files]):
        df,metrics,options = mt.ut.format_tuning(folder)
        L.append(df)
df = pd.concat(L)
df.loc[df['pp_method'] == 'mito_preprocessing', 'pp_method'] = 'MiTo'
df['combo'] = df[['pp_method', 'filtering']].astype(str).agg('-'.join, axis=1)

# Define options and metrics
varying_options = (df[options].nunique()).loc[lambda x:x>1].index.to_list()
metrics_of_interest = ['ARI', 'NMI', 'corr', 'AUPRC', 'n_cells', 'n_vars', 'n_GBC_groups']
metrics_of_interest += [
    'freq_lineage_biased_muts', 'median_n_vars_per_cell', 
    'transitions_vs_transversions_ratio', 'n_dbSNP', 'n_REDIdb'
]


##


# 1. Extended summary -------------------------- 

# Params
plu.set_rcParams()
matplotlib.rcParams.update({'figure.dpi':150})

# (
#     df.groupby(['sample', 'combo'])[metrics_of_interest]
#     .median()
# )


##


# PP method - filtering combo
fig =  plt.figure(figsize=(14,6.5))

x_order = ['MDA_clones', 'MDA_lung', 'MDA_PT']
by_order = ['samtools-baseline', 'freebayes-baseline', 'cellsnp-lite-MQuad', 'MiTo-MiTo', 'maegatk-MiTo']
cmap = plu.create_palette(df, 'combo', order=by_order, col_list=sc.pl.palettes.vega_10_scanpy)
test = (df['combo'].isin(['maegatk-MQuad', 'MiTo-MQuad']))

for i,metric in enumerate(metrics_of_interest):
    ax = fig.add_subplot(2,6,i+1)
    plu.bar(
        df.loc[~test].set_index('job_id'), 
        x='sample', y=metric, 
        by='combo',
        x_order=x_order,
        by_order=by_order,
        categorical_cmap=cmap,
        ax=ax
    )
    plu.format_ax(ax=ax, 
                  xticks=x_order if i>=6 else [ '' for _ in range(3) ], 
                  rotx=90,
                  xlabel='', ylabel=metric, reduced_spines=True)
    if i==5:
        plu.add_legend(cmap, label='Preprocessing-\nMT-SNVs filtering\ncombination', 
                       ax=ax, bbox_to_anchor=(1,0.3), loc='upper left')

fig.subplots_adjust(top=.90, bottom=.25, right=.82, left=.05, wspace=.8)
plt.show()


##


# 2. maegatk/MiTo with and without MiTo filter -------------------------- 

x_order = ['MDA_clones', 'MDA_lung', 'MDA_PT']
by_order = ['MiTo-MQuad', 'MiTo-MiTo', 'maegatk-MQuad', 'maegatk-MiTo']
colors = ['#EE8383', '#d62728', '#D99BEC', '#aa40fc']
cmap = dict(zip(by_order, colors))


##

fig, axs = plt.subplots(1,2,figsize=(7,3.5))

plu.bar(
    df.loc[df['combo'].isin(by_order)].set_index('job_id'), 
    x='sample', y='ARI', 
    by='combo',
    x_order=x_order,
    by_order=by_order,
    categorical_cmap=cmap,
    ax=axs[0]
)
plu.format_ax(ax=axs[0], xlabel='', ylabel='ARI', reduced_spines=True, rotx=90)
plu.bar(
    df.loc[df['combo'].isin(by_order)].set_index('job_id'), 
    x='sample', y='NMI', 
    by='combo',
    x_order=x_order,
    by_order=by_order,
    categorical_cmap=cmap,
    ax=axs[1]
)
plu.format_ax(ax=axs[1], xlabel='', ylabel='NMI', reduced_spines=True, rotx=90)
plu.add_legend(cmap, label='Preprocessing-\nMT-SNVs filtering\ncombination', ax=axs[1])

fig.subplots_adjust(top=.90, bottom=.3, right=.7, left=.1, wspace=.35)
plt.show()


##


# 3. MiTo-vs-maegatk  -------------------------- 

path_data = os.path.join(path_main, 'data', 'bench', 'MiTo_vs_maegatk_consensus')
path_QC_cells = os.path.join(path_main, 'data', 'general', 'QC', 'MDA_clones.txt')

##


def get_metrics(path):

    L = []
    for base in ['A', 'C', 'T', 'G']:

        print(f'base {base}...')
        base_df = pd.read_csv(os.path.join(path, f'{base}.txt.gz'), header=None)
        base_df.columns = ['pos', 'cell', 
                            'count_fw', 'qual_fw', 'cons_fw', 'gs_fw', 
                            'count_rev', 'qual_rev', 'cons_rev', 'gs_rev']

        mean_qual = np.nanmean(np.where(base_df[['qual_fw', 'qual_rev']]>0,base_df[['qual_fw', 'qual_rev']],np.nan), axis=1)
        mean_gs = np.nanmean(np.where(base_df[['gs_fw', 'gs_rev']]>0,base_df[['gs_fw', 'gs_rev']],np.nan), axis=1)
        mean_cons = np.nanmean(np.where(base_df[['cons_fw', 'cons_rev']]>0,base_df[['cons_fw', 'cons_rev']],np.nan), axis=1)

        df = pd.DataFrame({'pos':base_df['pos'], 'cell':base_df['cell'], 'qual':mean_qual, 'group_size':mean_gs, 'consensus':mean_cons })
        L.append(df)
    
    df = pd.concat(L)

    return df


##


# Gather coverage info
maegatk_mock = pd.read_csv(os.path.join(path_data, 'mock_maegatk', 'coverage.txt.gz'), header=None)
maegatk_mock.columns = ['pos', 'cell', 'n']
mito_prep = pd.read_csv(os.path.join(path_data, 'mito_preprocessing', 'coverage.txt.gz'), header=None)
mito_prep.columns = ['pos', 'cell', 'n']
# To wide
maegatk_mock_wide = maegatk_mock.pivot_table(index='cell', columns='pos', values='n').fillna(0)
mito_prep_wide = mito_prep.pivot_table(index='cell', columns='pos', values='n').fillna(0)
sites = list(set(mito_prep_wide.columns) & set(maegatk_mock_wide.columns))
maegatk_mock_wide = maegatk_mock_wide[sites]
mito_prep_wide = mito_prep_wide[sites]
# Gather qualified cells and subset both coverage table for them
cells = pd.read_csv(path_QC_cells, header=None)[0].unique().tolist()
cells = [ x.split('_')[0] for x in  cells ] 
subset = list(set(cells) & set(maegatk_mock_wide.index) & set(mito_prep_wide.index))
maegatk_mock_wide = maegatk_mock_wide.loc[subset]
mito_prep_wide = mito_prep_wide.loc[subset]

# Coverage metrics
cell_cov_mito_prep = mito_prep_wide.loc[:,mask_mt_sites(mito_prep_wide.columns)].median(axis=1)
cell_cov_maegatk = maegatk_mock_wide.loc[:,mask_mt_sites(mito_prep_wide.columns)].median(axis=1) 
frac_covered_sites_mito_prep = np.sum(mito_prep_wide.loc[:,mask_mt_sites(mito_prep_wide.columns)]>0, axis=1) / mask_mt_sites(mito_prep_wide.columns).sum()
frac_covered_sites_maegatk = np.sum(maegatk_mock_wide.loc[:,mask_mt_sites(maegatk_mock_wide.columns)]>0, axis=1) / mask_mt_sites(maegatk_mock_wide.columns).sum()


##


# Other metrics
df_mito_prep = get_metrics(os.path.join(path_data, 'mito_preprocessing')).query('cell in @subset')
df_maegatk = get_metrics(os.path.join(path_data, 'mock_maegatk')).query('cell in @subset')

# Basecalls
afm_mito_prep = sc.read(os.path.join(path_data, 'mito_preprocessing', 'afm.h5ad'))
afm_mito_prep = afm_mito_prep[[ f'{x}_MDA_clones' for x in subset],:].copy()
mito_prep_var_basecalls_per_cell = (afm_mito_prep.layers['AD'].A>0).sum(axis=1)
afm_maegatk = sc.read(os.path.join(path_data, 'mock_maegatk', 'afm.h5ad'))
afm_maegatk = afm_maegatk[[ f'{x}_MDA_clones' for x in subset],:].copy()
maegatk_var_basecalls_per_cell = (afm_maegatk.layers['AD'].A>0).sum(axis=1)


##


# Viz
fig, axs = plt.subplots(1,2,figsize=(10,5), subplot_kw={'projection': 'polar'})
mt.pl.MT_coverage_by_gene_polar(mito_prep.query('cell in @subset'), sample='MiTo', ax=axs[0])
mt.pl.MT_coverage_by_gene_polar(maegatk_mock.query('cell in @subset'), sample='MiTo', ax=axs[1])
fig.subplots_adjust(top=.8, bottom=.2, left=.2, right=.8 )
plt.show()



##


fig, axs = plt.subplots(1,6,figsize=(12,3))

# Median coverage per cell
ax = axs[0]
df = pd.concat([
    cell_cov_mito_prep.to_frame('cell_cov').assign(pp_method='MiTo'),
    cell_cov_maegatk.to_frame('cell_cov').assign(pp_method='maegatk'),
])
plu.box(df, x='pp_method', y='cell_cov', color='white', ax=ax, x_order=['maegatk','MiTo'])
plu.format_ax(ax=ax, reduced_spines=True, xlabel='', ylabel='Cell coverage')

# n covered sites per cell
ax = axs[1]
df = pd.concat([
    frac_covered_sites_mito_prep.to_frame('frac_sites').assign(pp_method='MiTo'),
    frac_covered_sites_maegatk.to_frame('frac_sites').assign(pp_method='maegatk'),
])
plu.box(df, x='pp_method', y='frac_sites', color='white', ax=ax, x_order=['maegatk','MiTo'])
plu.format_ax(ax=ax, reduced_spines=True, xlabel='', ylabel='Fraction of sites covered')

# Basecall UMI group size
ax = axs[2]
df = pd.concat([
    df_mito_prep.assign(pp_method='MiTo'),
    df_maegatk.assign(pp_method='maegatk')
])
plu.box(df, x='pp_method', y='group_size', color='white', ax=ax, x_order=['maegatk','MiTo'])
plu.format_ax(ax=ax, reduced_spines=True, xlabel='', ylabel='UMI group-size')

# Basecall quality
ax = axs[3]
plu.box(df, x='pp_method', y='qual', color='white', ax=ax, x_order=['maegatk','MiTo'])
plu.format_ax(ax=ax, reduced_spines=True, xlabel='', ylabel='Quality')
n_mito = np.sum(df_mito_prep['qual']<30)
n_maegatk = np.sum(df_maegatk['qual']<30)
ax.text(.1,.5, f'Qual<30', transform=ax.transAxes)
ax.text(.1,.43, f'MiTo: {n_mito}', transform=ax.transAxes)
ax.text(.1,.36, f'maegatk: {n_maegatk}', transform=ax.transAxes)

# Basecall consensus
ax = axs[4]
plu.box(df, x='pp_method', y='consensus', color='white', ax=ax, x_order=['maegatk','MiTo'])
plu.format_ax(ax=ax, reduced_spines=True, xlabel='', ylabel='Consensus score')
n_mito = np.sum(df_mito_prep['consensus']<.7)
n_maegatk = np.sum(df_maegatk['consensus']<.7)
ax.text(.1,.5, f'Consensus <0.7', transform=ax.transAxes)
ax.text(.1,.43, f'MiTo: {n_mito}', transform=ax.transAxes)
ax.text(.1,.36, f'maegatk: {n_maegatk}', transform=ax.transAxes)

# n variant basecalls per cell
ax = axs[5]
df = pd.concat([
    pd.Series(mito_prep_var_basecalls_per_cell).to_frame('n variant basecalls').assign(pp_method='MiTo'),
    pd.Series(maegatk_var_basecalls_per_cell).to_frame('n variant basecalls').assign(pp_method='maegatk'),
])
plu.box(df, x='pp_method', y='n variant basecalls', color='white', ax=ax, x_order=['maegatk','MiTo'])
plu.format_ax(ax=ax, reduced_spines=True, xlabel='', ylabel='n variant basecalls')

fig.tight_layout()
plt.show()

