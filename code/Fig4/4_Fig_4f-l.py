"""
Fig 4:
- Associate clonal features (size and branchedness) to gene expression
- Visualization of:
  1. Key clonal features (i.e., size and branchedness, Fig 4f)
  2. Association of clonal features to gene expression (Fig 4h-l)
"""

import os
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import mito as mt
import plotting_utils as plu
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_GEX = os.path.join(path_main, 'data', 'longitudinal')
path_lineage_inference = os.path.join(path_main, 'results', 'others', 'Fig4', 'lineage_inference')
path_fitness = os.path.join(path_main, 'results', 'others', 'Fig4', 'fitness')
path_figures = os.path.join(path_main, 'results', 'figures', 'Fig4')

# Set visualization params
plu.set_rcParams({'figure.dpi':350})


##


# 1. Fig 4f. Visualize clonal size and branchedness ------------------------------------#

# Read annotated cell phylogeny
with open(os.path.join(path_lineage_inference, 'annotated_tree.pickle'), 'rb') as f:
    tree = pickle.load(f)

# Retrieve clonal stats
df_cov = mt.ut.get_internal_node_stats(tree).loc[lambda x: x['clonal_node']]

# z-score fitness (LBI method by Neher et al., 2014, implemented in Cassiopeia)
df_cov['fitness'] = (df_cov['fitness']-df_cov['fitness'].mean()) / df_cov['fitness'].std()
df_cov = df_cov.sort_values('clade_size', ascending=False)

# Plot
fig, ax = plt.subplots(figsize=(3.5,1.5))
mt.pl.packed_circle_plot(
    df_cov, ax=ax, covariate='clade_size', cmap='Spectral_r', color_by='fitness'
)
covariate = df_cov['fitness']
plu.add_cbar(
    covariate, palette='Spectral_r', ax=ax, label='Branchedness', 
    layout=( (1.1,.5,.6,.06), 'top', 'horizontal' )
)
fig.subplots_adjust(left=.15, right=.5, top=.85, bottom=.1)
fig.savefig(os.path.join(path_figures, 'Fig_4f.pdf'))


##


# 2. NB regression to associate clonal features to gene expression --------------------#

# Read data
afm = sc.read(os.path.join(path_lineage_inference, 'afm_filtered.h5ad'))
adata = sc.read(os.path.join(path_GEX, 'expression.h5ad'))

# NB regression 
covariates = ['clade_size', 'fitness']

# Here we go:
for covariate in covariates:

    # NB regression: clone pseudobulk samples GEX ~ covariate
    agg = mt.tl.agg_pseudobulk(tree, adata, agg_method='sum', min_n_cells=30, n_cells=30)
    features = agg.columns[agg.columns.isin(adata.var_names)]
    results_nb = mt.tl.nb_regression(agg, features=features, predictor=covariate)
    results_nb.to_csv(os.path.join(path_fitness, f'nb_results_{covariate}.csv'))

    # GSEA
    results = results_nb.query('param==@covariate').set_index('gene')
    ranked_list = results['coef'].sort_values(ascending=False)
    gsea, gsea_df = mt.ut.run_GSEA(
        ranked_list, 
        collections=['GO_Biological_Process_2025', 'MSigDB_Hallmark_2020'],
        max_size_set=2000,
        min_size_set=15,
        max_pval_adj=.05
    )
    gsea_df['Term_adj'] = gsea_df['Term'].str.replace(
        r'(MSigDB_Hallmark_2020__|GO_Biological_Process_2025__|\(.*?\))',
        '',
        regex=True
    )
    gsea_df = gsea_df[['Term', 'Term_adj', 'ES', 'NES', 'pval_adj', 'Lead_genes']]
    gsea_df.to_csv(os.path.join(path_fitness, f'GSEA_{covariate}.csv'))

    # ORA
    gene_list = results['coef'].sort_values(ascending=False).head(50).index.to_list()
    ora, ora_df = mt.ut.run_ORA(
        gene_list, 
        collections=['GO_Biological_Process_2025', 'MSigDB_Hallmark_2020'],
        max_pval_adj=.05
    )
    ora_df['Term_adj'] = ora_df['Term'].str.replace(
        r'(MSigDB_Hallmark_2020__|GO_Biological_Process_2025__|\(.*?\))',
        '',
        regex=True
    )
    ora_df = ora_df[['Term', 'Term_adj', 'Overlap', 'Odds Ratio', 'pval_adj', 'Genes']]
    ora_df.to_csv(os.path.join(path_fitness, f'ORA_{covariate}.csv'))
 

##


# 2. Fig 4g-l. Visualization of NB regression results --------------------#

# NB regression results
fitness_nb = (
    pd.read_csv(os.path.join(path_fitness, f'nb_results_fitness.csv'), index_col=0)
    .set_index('gene')
    .drop(columns=['pval', '-logp10', 'param'])
    .rename(columns={'coef':'coef_fitness'})
)
clade_size_nb = (
    pd.read_csv(os.path.join(path_fitness, f'nb_results_clade_size.csv'), index_col=0)
    .set_index('gene')
    .drop(columns=['pval', '-logp10', 'param'])
    .rename(columns={'coef':'coef_clade_size'})
)
df = fitness_nb.join(clade_size_nb)
df['mean'] = df.mean(axis=1)


##


# Fig 4g
fig, ax = plt.subplots(figsize=(2.5,2.5))
plu.volcano(
    df, 
    'coef_fitness', 'coef_clade_size', 
    ax=ax, fig=fig, xlim=(-0.3,1), ylim=(-0.3,1), 
    cmap={'others':None, 'labelled':'r'},
    kwargs_labelled={'s':10},
    kwargs_text={'textsize':5, 'max_distance':.1}
)
ax.set(xlabel='Branchedness', ylabel='Clone size')
ax.set_xlim((-1.3,3))
ax.set_ylim((-1.3,3))
ax.axvline(0, color='r', linestyle='-', linewidth=.5)
ax.axhline(0, color='r', linestyle='-', linewidth=.5)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'Fig_4h.pdf'))


##


# GSEA mean fitted coefficients: key pathways significantly associated
ranked_list = df['mean'].sort_values(ascending=False)
gsea, gsea_df = mt.ut.run_GSEA(
    ranked_list, 
    collections=['GO_Biological_Process_2025', 'MSigDB_Hallmark_2020'],
    max_size_set=2000,
    min_size_set=15,
    max_pval_adj=.05
)
gsea_df['Term_adj'] = gsea_df['Term'].str.replace(
    r'(MSigDB_Hallmark_2020__|GO_Biological_Process_2025__|\(.*?\))',
    '',
    regex=True
)
gsea_df = gsea_df[['Term', 'Term_adj', 'ES', 'NES', 'pval_adj', 'Lead_genes']]
gsea_df.to_csv(os.path.join(path_fitness, f'GSEA_mean_coeffs.csv'))


##


# Fig 4h
fig, axs = plt.subplots(3,1,figsize=(2,3), sharex=True)

d_terms = {
    'MT-metabolism' : 'GO_Biological_Process_2025__Mitochondrial Gene Expression (GO:0140053)',
    'FAs oxidation' : 'GO_Biological_Process_2025__Fatty Acid Oxidation (GO:0019395)',
    'PDGF signaling' : 'GO_Biological_Process_2025__Platelet-Derived Growth Factor Receptor Signaling Pathway (GO:0048008)'
}

for i,name in enumerate(d_terms):
    term = d_terms[name]
    axs[i].plot(gsea.results[term]['RES'])
    fdr = gsea.results[term]['fdr']
    plu.format_ax(
        ax=axs[i], 
        title=name, 
        reduced_spines=True,
        ylabel='' if i!=1 else 'ES', 
        xlabel='' if i!=2 else 'Rank'
    )
    axs[i].text(.6, .05, f'FDR<{fdr:.1e}', transform=axs[i].transAxes, fontsize=6)

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'Fig_4h.pdf'))


##


# Selected 6-genes signature
adata = adata[tree.cell_meta.index].copy()
adata.obs['MiTo clone'] = tree.cell_meta['MiTo clone']
adata.obs['top_clones'] = np.where(
    adata.obs['MiTo clone'].isin(['MT-0', 'MT-1', 'MT-2']), 
    'top', 'others'
)
genes = ['C5orf46', 'CSAG1', 'DMKN', 'LSR', 'MAGEA12', 'MAGEA3']
sc.tl.score_genes(adata, gene_list=genes, score_name='6_genes')


##


# Fig 4i
df_genes = (
    pd.DataFrame(
        adata[:,genes].layers['raw'].toarray(), # Raw values, as in sc.pl.dotplot
        columns=genes, index=adata.obs_names
    )
    .join(adata.obs[['top_clones']])
)
df_plot = (   
    df_genes.groupby('top_clones').mean()
    .reset_index().melt(id_vars='top_clones', value_name='mean_expr', value_vars=genes)
    .merge(
        df_genes.groupby('top_clones').apply(lambda x: (x>0).sum()/x.shape[0])
        .reset_index().melt(id_vars='top_clones', value_name='freq_positive', value_vars=genes),
        on=['variable', 'top_clones']
    )
)

fig, ax = plt.subplots(figsize=(3.5,2))
plu.dotplot(df_plot, x='variable', y='top_clones', order_y=['top', 'others'],
            color='mean_expr', size='freq_positive', palette='viridis', ax=ax)
ax.margins(x=0.1, y=0.5)
plu.format_ax(ax=ax, xlabel='', ylabel='', rotx=90)
ax.get_legend().set_visible(False)
plu.add_cbar(x=df_plot['mean_expr'], palette='viridis', ax=ax, label='Gene expression', 
             layout=( (-0.2,1.5,.6,.2), 'top', 'horizontal' ))
fig.subplots_adjust(left=.35, right=.65, top=.6, bottom=.45)
fig.savefig(os.path.join(path_figures, 'Fig_4i.pdf'))


##


# Fig 4l
fig, axs = plt.subplots(1,2,figsize=(3,1.5))

cmap = {'top':'#E1899E', 'others':'#72C1AF'}
sc.pl.umap(
    adata, color='top_clones', ax=axs[0], 
    palette=cmap,
    show=False, legend_loc=None, 
    frameon=False, 
    title='Top clones',
    size=10,
)
plu.add_legend(
    ax=axs[0], colors=cmap, label='', 
    bbox_to_anchor=(1,0), loc='lower right', 
    artists_size=5, ticks_size=5
)
sc.pl.umap(
    adata, color='6_genes', ax=axs[1], 
    cmap='viridis',
    show=False, colorbar_loc=None, 
    frameon=False, 
    title='6-genes',
    size=10,
)
plu.add_cbar(
    x=adata.obs['6_genes'], 
    palette='viridis',
    ax=axs[1],
    label_size=5, 
    ticks_size=5,
    layout=( (.95,.05,.02,.22), 'left','vertical' )
)

fig.subplots_adjust(top=.85, bottom=.1, left=.1, right=.9, wspace=.005)
fig.savefig(os.path.join(path_figures, 'Fig_4l.pdf'))


##

