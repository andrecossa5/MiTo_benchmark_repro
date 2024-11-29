"""
Quantify differences ans similarity between maegatk and mito_preprocessing.
"""

import os
import pandas as pd
import numpy as np
from mito_utils.preprocessing import *
from mito_utils.plotting_base import *
matplotlib.use('macOSX')
set_rcParams()


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


# Paths
path_QC_cells = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/data/MI_TO_bench/QC_cells/MDA_clones.txt'
path_results = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/results/MI_TO_bench/maegatk_vs_mito'

# Gather coverage info
maegatk_mock = pd.read_csv(os.path.join(path_results, 'mock_maegatk', 'coverage.txt.gz'), header=None)
maegatk_mock.columns = ['pos', 'cell', 'n']
mito_prep = pd.read_csv(os.path.join(path_results, 'mito_preprocessing', 'coverage.txt.gz'), header=None)
mito_prep.columns = ['pos', 'cell', 'n']
# Wide
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

# np.corrcoef(maegatk_mock_wide.values.flatten(), mito_prep_wide.values.flatten())[0,1]

cell_cov_mito_prep = mito_prep_wide.loc[:,mask_mt_sites(mito_prep_wide.columns)].median(axis=1)
cell_cov_maegatk = maegatk_mock_wide.loc[:,mask_mt_sites(mito_prep_wide.columns)].median(axis=1) 
frac_covered_sites_mito_prep = np.sum(mito_prep_wide.loc[:,mask_mt_sites(mito_prep_wide.columns)]>0, axis=1) / mask_mt_sites(mito_prep_wide.columns).sum()
frac_covered_sites_maegatk = np.sum(maegatk_mock_wide.loc[:,mask_mt_sites(maegatk_mock_wide.columns)]>0, axis=1) / mask_mt_sites(maegatk_mock_wide.columns).sum()


##


# Other metrics
df_mito_prep = get_metrics(os.path.join(path_results, 'mito_preprocessing')).query('cell in @subset')
df_maegatk = get_metrics(os.path.join(path_results, 'mock_maegatk')).query('cell in @subset')

# np.sum(df_mito_prep['consensus']<.7) / df_mito_prep.shape[0]
# np.sum(df_maegatk['consensus']<.7) / df_maegatk.shape[0]
# np.sum(df_mito_prep['qual']<30) / df_mito_prep.shape[0]
# np.sum(df_maegatk['qual']<30) / df_maegatk.shape[0]

# Basecalls
afm_mito_prep = sc.read(os.path.join(path_results, 'mito_preprocessing', 'afm.h5ad'))
afm_mito_prep = afm_mito_prep[[ f'{x}_MDA_clones' for x in subset],:].copy()
mito_prep_var_basecalls_per_cell = (afm_mito_prep.layers['AD'].A>0).sum(axis=1)
afm_maegatk = sc.read(os.path.join(path_results, 'mock_maegatk', 'afm.h5ad'))
afm_maegatk = afm_maegatk[[ f'{x}_MDA_clones' for x in subset],:].copy()
maegatk_var_basecalls_per_cell = (afm_maegatk.layers['AD'].A>0).sum(axis=1)


##


# Viz
colors = {'mito_preprocessing':'#09C18A', 'maegatk':'#FF5500'}

fig, axs = plt.subplots(2,3,figsize=(9,6), sharex=True)

# Median coverage per cell
ax = axs[0,0]
df = pd.concat([
    cell_cov_mito_prep.to_frame('cell_cov').assign(pp_method='mito_preprocessing'),
    cell_cov_maegatk.to_frame('cell_cov').assign(pp_method='maegatk'),
])
box(df, x='pp_method', y='cell_cov', c=colors, ax=ax, order=['maegatk','mito_preprocessing'])
format_ax(ax=ax, reduced_spines=True, ylabel='Cell coverage', xlabel_size=10, ylabel_size=10, xticks_size=10)

# n covered sites per cell
ax = axs[0,1]
df = pd.concat([
    frac_covered_sites_mito_prep.to_frame('frac_sites').assign(pp_method='mito_preprocessing'),
    frac_covered_sites_maegatk.to_frame('frac_sites').assign(pp_method='maegatk'),
])
box(df, x='pp_method', y='frac_sites', c=colors, ax=ax, order=['maegatk','mito_preprocessing'])
format_ax(ax=ax, reduced_spines=True, ylabel='Fraction of sites covered', xlabel_size=10, ylabel_size=10, xticks_size=10)

# Basecall UMI group size
ax = axs[0,2]
df = pd.concat([
    df_mito_prep.assign(pp_method='mito_preprocessing'),
    df_maegatk.assign(pp_method='maegatk')
])
box(df, x='pp_method', y='group_size', c=colors, ax=ax, order=['maegatk','mito_preprocessing'])
format_ax(ax=ax, reduced_spines=True, ylabel='UMI group-size', xlabel_size=10, ylabel_size=10, xticks_size=10)

# Basecall quality
ax = axs[1,0]
box(df, x='pp_method', y='qual', c=colors, ax=ax, order=['maegatk','mito_preprocessing'])
format_ax(ax=ax, reduced_spines=True, ylabel='Quality', xlabel_size=10, ylabel_size=10, xticks_size=10)
n_mito = np.sum(df_mito_prep['qual']<30)
n_maegatk = np.sum(df_maegatk['qual']<30)
ax.text(.25,.5, f'n qual<30:', transform=ax.transAxes)
ax.text(.25,.43, f'-mito_preprocessing: {n_mito}', transform=ax.transAxes)
ax.text(.25,.36, f'-maegatk: {n_maegatk}', transform=ax.transAxes)

# Basecall consensus
ax = axs[1,1]
box(df, x='pp_method', y='consensus', c=colors, ax=ax, order=['maegatk','mito_preprocessing'])
format_ax(ax=ax, reduced_spines=True, ylabel='Consensus score', xlabel_size=10, ylabel_size=10, xticks_size=10)
n_mito = np.sum(df_mito_prep['consensus']<.7)
n_maegatk = np.sum(df_maegatk['consensus']<.7)
ax.text(.25,.5, f'n consensus <0.7:', transform=ax.transAxes)
ax.text(.25,.43, f'-mito_preprocessing: {n_mito}', transform=ax.transAxes)
ax.text(.25,.36, f'-maegatk: {n_maegatk}', transform=ax.transAxes)

# n variant basecalls per cell
ax = axs[1,2]
df = pd.concat([
    pd.Series(mito_prep_var_basecalls_per_cell).to_frame('n variant basecalls').assign(pp_method='mito_preprocessing'),
    pd.Series(maegatk_var_basecalls_per_cell).to_frame('n variant basecalls').assign(pp_method='maegatk'),
])
box(df, x='pp_method', y='n variant basecalls', c=colors, ax=ax, order=['maegatk','mito_preprocessing'])
format_ax(ax=ax, reduced_spines=True, ylabel='n variant basecalls', xlabel_size=10, ylabel_size=10, xticks_size=10)

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'maegatk_vs_mito.pdf'))


##