"""
Link clonal fitness with gene expression programs.
"""

import os
import pickle
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
path_others = os.path.join(path_main, 'results', 'others', 'Fig4')
path_figures = os.path.join(path_main, 'results', 'figures', 'Fig4')

# Read data
afm = sc.read(os.path.join(path_mt, 'afm_filtered.h5ad'))
tree_metrics = pd.read_csv(os.path.join(path_mt, 'tree_metrics.csv'), index_col=0)
with open(os.path.join(path_mt, 'annotated_tree.pickle'), 'rb') as f:
    tree = pickle.load(f)
adata = sc.read(os.path.join(path_expression, 'expression.h5ad'))


##


# Prep pseudobulk data table
agg = (
    pd.DataFrame(adata.layers['raw'].toarray(), 
                 index=adata.obs_names, columns=adata.var_names)
    .loc[tree.leaves]
    .assign(clone=tree.cell_meta['MiTo clone'])
    .groupby('clone')
    .sum()
)
# Add covariates, re-scaled
clone_fitness = (
    mt.tl.get_internal_node_stats(tree)
    .loc[lambda x: x['clonal_node']].reset_index(names='lca')
    .merge(tree.cell_meta[['MiTo clone', 'lca']], on='lca')
    .drop_duplicates()
    .set_index('MiTo clone')
    ['fitness']
)
agg['fitness'] = clone_fitness.loc[agg.index]
agg['fitness'] = (agg['fitness'] - agg['fitness'].mean()) / agg['fitness'].std()
agg['counts'] = agg.sum(axis=1) 
agg['counts'] = (agg['counts'] - agg['counts'].mean()) / agg['counts'].std()

# Negative Binomial regression
results = mt.tl.nb_regression(agg, features=adata.var_names, predictor='fitness + counts')


# GSEA
# res = mt.ut.run_GSEA(results['coef'].sort_values(ascending=False), collection='GO_Biological_Process_2025')
# res['Term'] = res['Term'].map(lambda x: x.replace('GO_Biological_Process_2025__', ''))

# import gseapy
# names = pd.Series(gseapy.get_library_name())
# names[names.str.contains('GO')]


##


# Viz NB regression results

# 1. Volcano plot
fig, ax = plt.subplots(figsize=(3,3))
plu.volcano(
    results.query('param=="fitness"').set_index('gene'), 
    'coef', '-logp10', 
    ax=ax, fig=fig, xlim=(-2,3.5), ylim=3, 
    cmap={'others':None, 'labelled':'r'}
)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'volcano_fitness_nbreg.pdf'))


##


# 2. ...

