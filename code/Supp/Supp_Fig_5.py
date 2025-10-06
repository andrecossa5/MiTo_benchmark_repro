"""
Supp Fig 5
MiTo (raw sequencing reads pre-processing) vs maegatk (full stat comparison).
"""

import os
import logging
import pandas as pd
import numpy as np
import scanpy as sc
import mito as mt
import matplotlib
import matplotlib.pyplot as plt
import plotting_utils as plu
from mito.ut import mask_mt_sites
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'data', 'bench', 'MiTo_vs_maegatk_consensus')
path_QC_cells = os.path.join(path_main, 'data', 'general', 'QC', 'MDA_clones.txt')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')


# Set visualization params
plu.set_rcParams({'figure.dpi':350})


##


# Utils
def get_metrics(path):

    L = []
    for base in ['A', 'C', 'T', 'G']:

        logging.info(f'base {base}...')
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


# 1. Extract all metrics ------------------------------------#

# Gather coverage info
maegatk_mock = pd.read_csv(os.path.join(path_data, 'mock_maegatk', 'coverage.txt.gz'), header=None)
maegatk_mock.columns = ['pos', 'cell', 'n']
mito_prep = pd.read_csv(os.path.join(path_data, 'mito_preprocessing', 'coverage.txt.gz'), header=None)
mito_prep.columns = ['pos', 'cell', 'n']
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
frac_covered_sites_mito_prep = np.sum(
    mito_prep_wide.loc[:,mask_mt_sites(mito_prep_wide.columns)]>0, axis=1) / \
    mask_mt_sites(mito_prep_wide.columns).sum()
frac_covered_sites_maegatk = np.sum(
    maegatk_mock_wide.loc[:,mask_mt_sites(maegatk_mock_wide.columns)]>0, axis=1) / \
    mask_mt_sites(maegatk_mock_wide.columns).sum()


##


# Other metrics
df_mito_prep = get_metrics(os.path.join(path_data, 'mito_preprocessing')).query('cell in @subset')
df_maegatk = get_metrics(os.path.join(path_data, 'mock_maegatk')).query('cell in @subset')

# Basecall metrics

# mito_preprocessing
afm_mito_prep = mt.io.make_afm(
    os.path.join(path_data, 'mito_preprocessing'), 
    sample='', 
    pp_method='mito_preprocessing'
)
afm_mito_prep.obs_names = afm_mito_prep.obs_names.str.split('_').str[0] # Correct automatic renaming
afm_mito_prep = afm_mito_prep[subset,:].copy()
mito_prep_var_basecalls_per_cell = (afm_mito_prep.layers['AD'].toarray()>0).sum(axis=1)

# mock_maegatk
afm_maegatk = mt.io.make_afm(
    os.path.join(path_data, 'mock_maegatk'), 
    sample='', 
    pp_method='mito_preprocessing'
)
afm_maegatk.obs_names = afm_maegatk.obs_names.str.split('_').str[0]     # Correct automatic renaming
afm_maegatk = afm_maegatk[subset,:].copy()
maegatk_var_basecalls_per_cell = (afm_maegatk.layers['AD'].toarray()>0).sum(axis=1)


##


# 1. Visualize all metrics ------------------------------------#

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
fig.savefig(os.path.join(path_figures, 'Supp_Fig_5.pdf'))


##