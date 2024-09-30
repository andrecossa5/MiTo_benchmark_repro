"""
Additional quality check on mito_preprocessing output. MT- library.
"""

import os
import numpy as np
import pandas as pd
from mito_utils.preprocessing import *
from mito_utils.plotting_base import *
from mito_utils.diagnostic_plots import *
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
afm = make_AFM(path_afm, path_meta=None, sample=sample, nmads=5, mean_cov_all=25)

# Radial plot
fig, ax = plt.subplots(figsize=(4.5,4.5))
plot_ncells_nAD(afm, ax=ax, title=sample, s=5)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'MT_QC', f'{sample}_ncells_nAD.png'), dpi=1000)


fig, ax = plt.subplots(figsize=(4.5,4.5), subplot_kw={'projection': 'polar'})
MT_coverage_by_gene_polar(afm, ax=ax, sample=sample)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'MT_QC', f'{sample}_MT_coverage.png'), dpi=1000)


##


# Legend
df_mt = pd.DataFrame(MAESTER_genes_positions, columns=['gene', 'start', 'end']).set_index('gene').sort_values('start')
colors = { k:v for k,v in zip(df_mt.index, sc.pl.palettes.default_102[:df_mt.shape[0]])}

fig, ax = plt.subplots(figsize=(10,5))
ax.axis('off')
add_legend(ax=ax, label='Target MT-genes', colors=colors, ncols=round(len(colors)/2), bbox_to_anchor=(.5,.5), loc='center',
           label_size=12, artists_size=10, ticks_size=10)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'MT_QC', 'MT_genes.png'), dpi=1000)


##


# Expression of MT-genes
sample = 'MDA_clones'
path_afm = os.path.join(path_data, 'AFMs', 'MDA_clones.h5ad')
path_meta = os.path.join(path_data, 'cells_meta.csv')
afm = make_AFM(path_afm, path_meta=None, sample=sample, nmads=5, mean_cov_all=25)
mt_expr = pd.read_csv(os.path.join(path_data, 'miscellanea', 'mt_genes_expr.csv'), index_col=0)
cells = list( set(afm.obs_names) & set(mt_expr.index))

# MT-gene mean expression vs MAESTER mean base coverage
mean_expr = mt_expr.loc[cells,:].mean(axis=0)

# Annotate MT-genome sites
mt_genes_positions = [ x for x in all_mt_genes_positions if x[0] in mt_expr.columns ]
sites = afm.uns['per_position_coverage'].columns
annot = {}
for x in sites:
    x = int(x)
    mapped = False
    for mt_gene, start, end in mt_genes_positions:
        if x>=start and x<=end:
            annot[str(x)] = mt_gene
            mapped = True
    if not mapped:
        annot[str(x)] = 'other'

mean_site_cov = afm[cells,:].uns['per_position_coverage'].T.mean(axis=1).to_frame('cov')
mean_site_cov['gene'] = pd.Series(annot)
mean_site_cov = mean_site_cov.query('gene!="other"').groupby('gene')['cov'].mean()
mean_site_cov = mean_site_cov[mean_expr.index]

##

fig, ax = plt.subplots(figsize=(4.5,4.5))
ax.plot(mean_expr.values, mean_site_cov.values, 'ko')
sns.regplot(data=pd.DataFrame({'expr':mean_expr, 'cov':mean_site_cov}), 
            x='expr', y='cov', ax=ax, scatter=False)
format_ax(ax, title='MT-transcripts counts vs site coverage',
          xlabel='Mean expression (gene nUMI, 10x)', 
          ylabel='Mean site coverage (per-site nUMI, MAESTER)')
corr = np.corrcoef(mean_expr.values, mean_site_cov.values)[0,1]
ax.text(.05, .9, f'Pearson\'s r: {corr:.2f}', transform=ax.transAxes)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'MT_QC', f'{sample}_site_coverage_by_gene_expression.png'), dpi=500)


##