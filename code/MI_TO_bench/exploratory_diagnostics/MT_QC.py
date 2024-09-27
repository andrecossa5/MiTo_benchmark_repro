"""
Additional quality check on mito_preprocessing output. MT- library.
"""

import os
import numpy as np
import pandas as pd
from mito_utils.preprocessing import *
from mito_utils.plotting_base import *
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro'
path_data = os.path.join(path_main, 'data', 'MI_TO_bench')
path_results = os.path.join(path_main, 'results',  'MI_TO_bench')


##


# Plot individual metrics distribution using /tables of one sample 
feat = 'cons'

L = []
for base in ['A','C','T','G']:
    df = pd.read_csv(os.path.join(path_data, 'MT_tables', f'{base}.txt.gz'), header=None)
    df.columns = ['pos', 'cell', 'n_fw', 'q_fw', 'cons_fw', 'gs_fw', 'n_rev', 'q_rev', 'cons_rev', 'gs_rev'] 
    x = df[f'{feat}_fw'].loc[lambda x: x>0].to_list() + df[f'{feat}_rev'].loc[lambda x: x>0].to_list()
    L += x

x = pd.Series(np.array(L))

fig, ax = plt.subplots(figsize=(4,4))
sns.kdeplot(x, ax=ax, fill=True, alpha=.8)
format_ax(ax=ax, xlabel='Consensus UMI consensus score', ylabel='Density', reduced_spines=True)
ax.text(.1, .9, f'Median: {np.median(x)}', transform=ax.transAxes)
ax.text(.1, .85, f'Min-max: {np.min(x)}-{np.max(x)}', transform=ax.transAxes)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'MT_QC', f'{feat}_distribution.png'))


##


FW = []
REV = []
N = []
for base in ['A','C','T','G']:
    df = pd.read_csv(os.path.join(path_data, 'MT_tables', f'{base}.txt.gz'), header=None)
    df.columns = ['pos', 'cell', 'n_fw', 'q_fw', 'cons_fw', 'gs_fw', 'n_rev', 'q_rev', 'cons_rev', 'gs_rev'] 
    FW += df['n_fw'].to_list()
    REV += df['n_rev'].to_list()
    n = df.groupby('cell').apply(lambda x: np.sum((x['n_fw']>0) & (x['n_rev']>0)) / x['pos'].nunique()).to_list()
    N += n

FW = pd.Series(np.array(FW))
REV = pd.Series(np.array(REV))
N = pd.Series(np.array(N))


fig, ax = plt.subplots(figsize=(4,4))
ax.plot(FW, REV, 'o', alpha=.2)
format_ax(ax=ax, xlabel='n UMIs forward strand', ylabel='n UMIs reverse strand', reduced_spines=True)
n = N.median()
ax.text(.2, .85, f'Median (across cell) % sites\nwith bi-allelic expression: {n*100:.2f}%', transform=ax.transAxes)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'MT_QC', 'strand_bias.png'))


##


# Cell and site coverage
sample = 'MDA_clones'
path_afm = os.path.join(path_data, 'AFMs', 'MDA_clones.h5ad')
afm = read_one_sample(path_afm, path_meta=None, sample=sample, nmads=5, mean_coverage=25)

compute_metrics_raw(afm)


##


