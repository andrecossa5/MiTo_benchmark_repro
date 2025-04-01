"""
Main Fig.2. Tuning MT-SNVs workflow hyper-parameters, example output (MDA_PT)

1. Prepare grouped dataframes for funky heatmaps (all samples: MDA_clones, MDA_PT, MDA_lung).
2. Choose final MT-SNVs space for clonal inference benchmark.
3. Visualize example output: MDA_PT 
"""

import os
import pickle
import pandas as pd
import scanpy as sc
import mito as mt
import plotting_utils as plu
from itertools import chain
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_figures = os.path.join(path_main, 'results', 'figures', 'Fig2')
path_results = os.path.join(path_main, 'results', 'others', 'Fig2')


##


# 1. Group jobs and score --------------------------------------- # 

# Format tune output
path_data = os.path.join(path_main, 'data', 'tune')
df, metrics, options = mt.ut.format_tuning(path_data)

# Set annot
groupings = ['pp_method', 'af_confident_detection', 'min_n_confidently_detected', 'bin_method', 'metric']
metric_annot = {
    'Mutation Quality' : ['n_dbSNP', 'n_REDIdb', 'transitions_vs_transversions_ratio'],
    'Association with GBC' : ['freq_lineage_biased_muts', 'AUPRC', 'ARI', 'NMI'],                               
    'Tree structure' : ['corr'],
    'Connectedness' : ['density', 'transitivity', 'average_path_length', 'average_degree', 'proportion_largest_component'],
    'Variation' : ['genomes_redundancy', 'median_n_vars_per_cell'],                                                           
    'Yield' : ['n_GBC_groups', 'n_cells', 'n_vars']                                                                
}  
relevant_metrics = list(chain.from_iterable([ metric_annot[k] for k in metric_annot ]))
relevant_metrics = [ f'{x}_rescaled' for x in relevant_metrics ]
weights = {
    'Mutation Quality': .1,
    'Association with GBC': .4,
    'Tree structure' : .1,
    'Connectedness' : .0,
    'Variation' : 0,
    'Yield' : .4
}

# Score and rank, for each single task (i.e., sample)
n = 5
samples = ['MDA_clones', 'MDA_lung', 'MDA_PT']
for sample in samples:
    df_ = df.query('sample==@sample')
    df_ranked = mt.ut.rank_items(df_, groupings, metrics, weights, metric_annot)
    df_final = pd.concat([df_ranked.head(n), df_ranked.tail(n)])
    metric_type_scores = df_final.columns[df_final.columns.str.contains('score')].to_list()
    (
        df_final[groupings+metric_type_scores+relevant_metrics]
        .to_csv(os.path.join(path_results, f'{sample}_grouped.csv'))
    )

##


# 2. Choose final MT-SNVs spaces for clonal inference benchmark --------------------------------------- # 

# Format tune output
path_data = os.path.join(path_main, 'data', 'tune')
df, metrics, options = mt.ut.format_tuning(path_data)
df = df.query('bin_method!="MiTo_smooth"')

# Choose jobs: MDA_clones
sample = 'MDA_clones'
df_selected = (
    df.query('sample==@sample and n_cells>=250 and n_GBC_groups>=6 and n_vars>10')
)
df_selected.to_csv(os.path.join(path_results, f'{sample}_filtered_jobs.csv'))
df_selected = (
    df_selected[[
        'job_id', 'pp_method', 'bin_method', 'af_confident_detection', 'min_cell_number', 'metric',
        'ARI', 'corr', 'NMI', 'AUPRC', 'n_cells', 'unassigned', 'n_vars', 'n_GBC_groups', 'n MiTo clone',
    ]]
)

# Choose appropriate bins
df_selected['n_vars'].describe()
bins = [0,20,30,df_selected['n_vars'].max()]
df_selected['cut'] = pd.cut(df_selected['n_vars'], bins=bins)
df_selected['cut'].value_counts()

# Select jobs
df_selected = df_selected.groupby('cut').apply(lambda x: x.sort_values('ARI', ascending=False).head(5))[[
    'job_id', 'pp_method', 'bin_method', 'metric', 'af_confident_detection', 'min_cell_number',
    'ARI', 'corr', 'NMI', 'AUPRC', 'n_cells', 'unassigned', 'n_vars', 'n_GBC_groups', 'n MiTo clone', 
]]

# Drop duplicates
df_selected_ = df_selected[['n_cells', 'n_vars', 'unassigned', 'n MiTo clone', 'n_GBC_groups']].drop_duplicates()
df_selected = df_selected.loc[df_selected_.index]
df_selected


##


# Choose jobs: MDA_PT
sample = 'MDA_PT'
df_selected = (
    df.query('sample==@sample and n_cells>=850 and n_GBC_groups>=30 and n_vars>10')
)
df_selected.to_csv(os.path.join(path_results, f'{sample}_filtered_jobs.csv'))
df_selected = (
    df_selected[[
        'job_id', 'pp_method', 'bin_method', 'af_confident_detection', 'min_cell_number', 'metric',
        'ARI', 'corr', 'NMI', 'AUPRC', 'n_cells', 'unassigned', 'n_vars', 'n_GBC_groups', 'n MiTo clone',
    ]]
)

# Choose appropriate bins
df_selected['n_vars'].describe()
bins = [0,45,75,df_selected['n_vars'].max()]
df_selected['cut'] = pd.cut(df_selected['n_vars'], bins=bins)
df_selected['cut'].value_counts()

# Select jobs
df_selected = df_selected.groupby('cut').apply(lambda x: x.sort_values('ARI', ascending=False).head(5))[[
    'job_id', 'pp_method', 'bin_method', 'metric', 'af_confident_detection', 'min_cell_number',
    'ARI', 'corr', 'NMI', 'AUPRC', 'n_cells', 'unassigned', 'n_vars', 'n_GBC_groups', 'n MiTo clone', 
]]

# Drop duplicates
df_selected_ = df_selected[['n_cells', 'n_vars', 'unassigned', 'n MiTo clone', 'n_GBC_groups']].drop_duplicates()
df_selected = df_selected.loc[df_selected_.index]
df_selected


##


# Choose jobs: MDA_lung
sample = 'MDA_lung'
df_selected = (
    df.query('sample==@sample and n_cells>=850 and n_GBC_groups>=10 and n_vars>10')
)
df_selected.to_csv(os.path.join(path_results, f'{sample}_filtered_jobs.csv'))
df_selected = (
    df_selected[[
        'job_id', 'pp_method', 'bin_method', 'af_confident_detection', 'min_cell_number', 'metric',
        'ARI', 'corr', 'NMI', 'AUPRC', 'n_cells', 'unassigned', 'n_vars', 'n_GBC_groups', 'n MiTo clone',
    ]]
)

# Choose appropriate bins
df_selected['n_vars'].describe()
bins = [0,20,30,df_selected['n_vars'].max()]
df_selected['cut'] = pd.cut(df_selected['n_vars'], bins=bins)
df_selected['cut'].value_counts()

# Select jobs
df_selected = df_selected.groupby('cut').apply(lambda x: x.sort_values('ARI', ascending=False).head(5))[[
    'job_id', 'pp_method', 'bin_method', 'metric', 'af_confident_detection', 'min_cell_number',
    'ARI', 'corr', 'NMI', 'AUPRC', 'n_cells', 'unassigned', 'n_vars', 'n_GBC_groups', 'n MiTo clone', 
]]

# Drop duplicates
df_selected_ = df_selected[['n_cells', 'n_vars', 'unassigned', 'n MiTo clone', 'n_GBC_groups']].drop_duplicates()
df_selected = df_selected.loc[df_selected_.index]
df_selected


##


# 3. Plot representative example: MDA_PT --------------------------------------- # 

# See also explore output. ./results/others/explore_final_jobs
path_data = os.path.join(path_main, 'data', 'bench', 'clonal_inference')

# Choose MT-SNVs space
sample = 'MDA_PT'
job_id = '74466a5d58'
path_data = os.path.join(path_data, sample, job_id)

# Load: afm, tree, metrics, phylocorr
afm = sc.read(os.path.join(path_data, 'afm.h5ad'))
with open(os.path.join(path_data, 'annotated_tree.pickle'), 'rb') as f:
    tree = pickle.load(f)
metrics = pd.read_csv(os.path.join(path_data, 'tree_metrics.csv'), index_col=0)

# Load colors
path_colors = os.path.join(path_main, 'data', 'general')
with open(os.path.join(path_colors, 'clones_colors_sc.pickle'), 'rb') as f:
    colors = pickle.load(f)


##


# Viz
plu.set_rcParams()

# UMAP
mt.pp.reduce_dimensions(afm)

# Fig 2, bottom
fig, axs = plt.subplots(1,4,figsize=(16,4))

mt.pl.draw_embedding(afm, feature='GBC', ax=axs[0], categorical_cmap=colors)
mt.pl.heatmap_variants(afm, tree=tree, ax=axs[1])
mt.pl.heatmap_distances(afm, tree=tree, ax=axs[2])
cmaps = { 
    'GBC' : colors,
    'MiTo clone' : plu.create_palette(tree.cell_meta, 'MiTo clone', col_list=sc.pl.palettes.default_102)
}
# tree.cell_meta.loc[tree.cell_meta['MiTo clone'].isna(), 'MiTo clone'] = 'unassigned'
clonal_nodes = tree.cell_meta['lca'].unique()[1:]
mt.pl.plot_tree(
    tree, 
    features=['GBC', 'MiTo clone'], 
    categorical_cmaps=cmaps, 
    ax=axs[3], 
    colorstrip_width=5,
    internal_node_subset=clonal_nodes,
    feature_internal_nodes='support',
    internal_node_kwargs={'markersize':7.5}
)
tree_stats = mt.tl.get_internal_node_stats(tree)
plu.add_cbar(
    tree_stats['support'], palette='Spectral_r', ax=axs[3], 
    vmin=.5, vmax=.8, label='Support', layout='outside'
)

fig.tight_layout()
fig.savefig(os.path.join(path_figures, f'{sample}_{job_id}_visualization.png'), dpi=1000)


##