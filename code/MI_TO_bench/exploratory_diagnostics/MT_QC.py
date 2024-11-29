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
set_rcParams()


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro'
path_data = os.path.join(path_main, 'data', 'MI_TO_bench')
path_results = os.path.join(path_main, 'results',  'MI_TO_bench')


##


# Biases
fig, axs = plt.subplots(1,2,figsize=(7,3.5))

# MDA_clones
path_afm = os.path.join(path_data, 'MT_tables', 'afm.h5ad')
afm = sc.read(path_afm)
afm = filter_cells(afm, cell_filter='filter2')

# Strand bias
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

ax = axs[0]
ax.plot(FW, REV, 'o', alpha=.2)
format_ax(ax=ax, xlabel='n UMIs forward strand', ylabel='n UMIs reverse strand', reduced_spines=True, xlabel_size=10, ylabel_size=10)
n = N.median()
ax.text(.2, .85, f'Median (across cell) % sites\nwith bi-allelic expression: {n*100:.2f}%', transform=ax.transAxes)

##

# Expression bias
mt_cov = pd.read_csv(os.path.join(path_data, 'MT_tables', 'coverage.txt.gz'), header=None)
mt_cov.columns = ['pos', 'cell', 'n']
mt_cov = mt_cov.pivot_table(index='cell', columns='pos', values='n').fillna(0)
mt_cov.index = mt_cov.index.map(lambda x: f'{x}_MDA_clones')
mt_expr = pd.read_csv(os.path.join(path_data, 'miscellanea', 'mt_genes_expr.csv'), index_col=0)
cells = list( set(mt_cov.index) & set(mt_expr.index) )
mt_cov = mt_cov.loc[cells,:]
mean_expr = mt_expr.loc[cells,:]

# MT-gene mean expression vs MAESTER mean base coverage
mean_expr = mt_expr.mean(axis=0)
mt_genes_positions = [ x for x in all_mt_genes_positions if x[0] in mt_expr.columns ]
sites = mt_cov.columns
annot = {}
for x in sites:
    x = int(x)
    mapped = False
    for mt_gene, start, end in mt_genes_positions:
        if x>=start and x<=end:
            annot[x] = mt_gene
            mapped = True
    if not mapped:
        annot[x] = 'other'

mean_site_cov = mt_cov.mean(axis=0).to_frame('cov')
mean_site_cov['gene'] = mean_site_cov.index.map(annot)
mean_site_cov = mean_site_cov.query('gene!="other"').groupby('gene')['cov'].mean()
mean_site_cov = mean_site_cov[mean_expr.index]

ax = axs[1]
ax.plot(mean_expr.values, mean_site_cov.values, 'ko')
sns.regplot(data=pd.DataFrame({'expr':mean_expr, 'cov':mean_site_cov}), x='expr', y='cov', ax=ax, scatter=False)
format_ax(ax, xlabel='Mean expression (gene nUMI, 10x)', ylabel='Mean site coverage\n(per-site nUMI, MAESTER)', reduced_spines=True, xlabel_size=10, ylabel_size=10)
corr = np.corrcoef(mean_expr.values, mean_site_cov.values)[0,1]
ax.text(.05, .9, f'Pearson\'s r: {corr:.2f}', transform=ax.transAxes)

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'MT_QC', 'strand_and_expression_bias.png'), dpi=500)


##


# Cell and site coverage
samples = ['MDA_clones', 'MDA_PT', 'MDA_lung']

fig, axs = plt.subplots(1,3, figsize=(12,4), subplot_kw={'projection': 'polar'})

for sample,ax in zip(samples, axs):
    cov = pd.read_csv(os.path.join(path_data, 'AFMs', 'mito_preprocessing', sample, 'coverage.txt.gz'), header=None)
    cov.columns = ['pos', 'cell', 'n']
    cov['cell'] = cov['cell'].map(lambda x: f'{x}_{sample}')
    cells = pd.read_csv(os.path.join(path_data, 'QC_cells', f'{sample}.txt'), header=None)[0]
    cov = cov.query('cell in @cells')
    MT_coverage_by_gene_polar(cov, sample=sample, subset=cells, ax=ax)

fig.subplots_adjust(left=.1, right=.9, top=.7, bottom=.1, wspace=.4)
fig.savefig(os.path.join(path_results, 'MT_QC', 'MT_site_coverage.pdf'))


##


# Legend
df_mt = pd.DataFrame(MAESTER_genes_positions, columns=['gene', 'start', 'end']).set_index('gene').sort_values('start')
colors = { k:v for k,v in zip(df_mt.index, sc.pl.palettes.default_102[:df_mt.shape[0]])}

fig, ax = plt.subplots(figsize=(10,5))
ax.axis('off')
add_legend(ax=ax, label='Target MT-genes', colors=colors, ncols=round(len(colors)/2), bbox_to_anchor=(.5,.5), loc='center',
           label_size=12, artists_size=10, ticks_size=10)

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'MT_QC', 'MT_genes.pdf'))


##


# Radial plot
fig, axs = plt.subplots(1,3,figsize=(10,3.7))

xticks = [0,1,2,5,10,20,40,80,200,500] 
for i,sample in enumerate(samples):
    ax =  axs.ravel()[i]
    afm = sc.read(os.path.join(path_data, 'AFMs', 'mito_preprocessing', sample, 'afm.h5ad'))
    plot_ncells_nAD(afm, ax=ax, title=sample, s=3, xticks=xticks)

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'MT_QC', f'unfiltered_ncells_AD.png'), dpi=500)


##