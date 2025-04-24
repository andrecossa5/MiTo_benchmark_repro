"""
Link clonal fitness with gene expression programs.
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
path_figures = os.path.join(path_main, 'results', 'figures', 'Fig4')


##


# Paths
path_expression = os.path.join(path_main, 'data', 'longitudinal')
path_mt = os.path.join(path_main, 'results', 'others', 'Fig4', 'longitudinal')
path_others = os.path.join(path_main, 'results', 'others', 'Fig4', 'fitness')
path_figures = os.path.join(path_main, 'results', 'figures', 'Fig4')

# Read data
afm = sc.read(os.path.join(path_mt, 'afm_filtered.h5ad'))
tree_metrics = pd.read_csv(os.path.join(path_mt, 'tree_metrics.csv'), index_col=0)
with open(os.path.join(path_mt, 'annotated_tree.pickle'), 'rb') as f:
    tree = pickle.load(f)
adata = sc.read(os.path.join(path_expression, 'expression.h5ad'))


##


# Viz clone branchedness and size
plu.set_rcParams()

df_cov = mt.ut.get_internal_node_stats(tree).loc[lambda x: x['clonal_node']]
df_cov['fitness'] = (df_cov['fitness']-df_cov['fitness'].mean()) / df_cov['fitness'].std()
df_cov = df_cov.sort_values('clade_size', ascending=False)

fig, ax = plt.subplots(figsize=(3.5,1.5))
mt.pl.packed_circle_plot(df_cov, ax=ax, covariate='clade_size', cmap='Spectral_r', color_by='fitness')
plu.add_cbar(df_cov['fitness'], palette='Spectral_r', ax=ax, label='Branchedness',
             layout=( (1.1,.5,.6,.06), 'top', 'horizontal' ))
fig.subplots_adjust(left=.15, right=.5, top=.85, bottom=.1)
fig.savefig(os.path.join(path_figures, 'bubble_clones.pdf'))


##


# Pseudobulk NB regression workflow
param = 'clade_size'
agg = mt.tl.agg_pseudobulk(tree, adata, agg_method='sum', min_n_cells=30, n_cells=30)
features = agg.columns[agg.columns.isin(adata.var_names)]
results_nb = mt.tl.nb_regression(agg, features=features, predictor=f'{param}')
results_nb.to_csv(os.path.join(path_others, f'nb_results_{param}.csv'))

# Gene ranks interpretations
results = results_nb.query('param==@param').set_index('gene')
    
# GSEA
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
gsea_df.to_csv(os.path.join(path_others, f'GSEA_{param}.csv'))

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
ora_df.to_csv(os.path.join(path_others, f'ORA_{param}.csv'))
 

##


# Visialization NB regression results
plu.set_rcParams({'figure.dpi':150})

# Read NB regression results
fitness_nb = (
    pd.read_csv(os.path.join(path_others, f'nb_results_fitness.csv'), index_col=0)
    .set_index('gene')
    .drop(columns=['pval', '-logp10', 'param'])
    .rename(columns={'coef':'coef_fitness'})
)
clade_size_nb = (
    pd.read_csv(os.path.join(path_others, f'nb_results_clade_size.csv'), index_col=0)
    .set_index('gene')
    .drop(columns=['pval', '-logp10', 'param'])
    .rename(columns={'coef':'coef_clade_size'})
)
df = fitness_nb.join(clade_size_nb)
df['mean'] = df.mean(axis=1)
# df['mean'].loc[lambda x: x>1.5].index
# df.sort_values('mean', ascending=False).head(20)


##


# 1. Volcano plot, fitted coefficients
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
fig.savefig(os.path.join(path_figures, 'volcano_fitness_nbreg_size.pdf'))


##


# 2. GSEA mean fitted coefficients
ranked_list = df.mean(axis=1).sort_values(ascending=False)
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
gsea_df.to_csv(os.path.join(path_others, f'GSEA_mean_coeffs.csv'))

# gsea_df[['Term', 'ES', 'NES', 'pval_adj']].head(10)['Term'].to_list()
# gsea_df[['Term_adj','NES']]
# gsea_df['Term'].to_list()

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
fig.savefig(os.path.join(path_figures, 'GSEA_fitness_nbreg.pdf'))


##


# 3. Selected 6-genes signature, dotplot and umaps
adata = adata[tree.cell_meta.index].copy()
adata.obs['MiTo clone'] = tree.cell_meta['MiTo clone']
adata.obs['top_clones'] = np.where(
    adata.obs['MiTo clone'].isin(['MT-0', 'MT-1', 'MT-2']), 
    'top', 'others'
)
genes = ['C5orf46', 'CSAG1', 'DMKN', 'LSR', 'MAGEA12', 'MAGEA3']
sc.tl.score_genes(adata, gene_list=genes, score_name='6_genes')

##

# Dotplot
df_genes = (
    pd.DataFrame(
        adata[:,genes].layers['raw'].toarray(), 
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
fig.savefig(os.path.join(path_figures, 'dotplot_genes_fitness.pdf'))

##

# UMAPs
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
fig.savefig(os.path.join(path_figures, 'umaps_gene_fitness.pdf'))


##

