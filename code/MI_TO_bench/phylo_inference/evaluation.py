"""
Evaluate phylogeny inference output.
"""

import os
from mito_utils.utils import *
from mito_utils.distances import *
from mito_utils.dimred import *
from mito_utils.clustering import *
from mito_utils.phylo import *
from mito_utils.phylo import *
from mito_utils.plotting_base import *
from mito_utils.embeddings_plots import *
from mito_utils.diagnostic_plots import *
from mito_utils.phylo_plots import *
matplotlib.use('macOSX')


##


# Get metrics
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro'
path_results = os.path.join(path_main, 'results', 'MI_TO_bench', 'phylo_inference')

# Retrieve and collapse everything
metrics = ['cell_filter', 'max_AD_counts', 'min_AD', 'min_mean_AD_in_confident', 'solver']    # The ones varied
          # 'min_median_af', 'min_n_positive', 'min_n_confidently_detected', 'distance_metric', 'af_confident_detection']

# Muts
df = get_metrics(path_results, 'dataset_df')
df = df.rename(columns={'metric_x':'metric', 'metric_y':'distance_metric'})[['metric', 'metric_value', 'sample']+metrics]
df.groupby(metrics+['metric'])['metric_value'].describe()
l = ['median_n_vars_per_cell', 'density', 'transitions_vs_transversions_ratio', 'freq_lineage_biased_muts', 'average_degree', 'transitivity', 'n_cells', 'n_vars']
df.query('metric in @l').groupby(['sample','metric'])['metric_value'].max()

# Distances
df = get_metrics(path_results, 'distance_metrics')
df = df.rename(columns={'metric_x':'metric', 'metric_y':'distance_metric'})[['metric', 'metric_value', 'sample']+metrics]
df.groupby(metrics+['metric'])['metric_value'].describe()
df.groupby(['sample','metric'])['metric_value'].max()

# Trees
df = get_metrics(path_results, 'tree_metrics')
df = df.rename(columns={'metric_x':'metric', 'metric_y':'distance_metric'})[['metric', 'metric_value', 'sample']+metrics]
df.groupby(metrics+['metric'])['metric_value'].describe()
l = ['median_support_upmost_nodes', 'corr_distances', 'total_assigned_char', 'median_n_assigned_char_per_clone', 'median_CI', 'median_RI', 'ARI', 'NMI']
df.query('metric in @l').groupby(['sample','metric'])['metric_value'].max()


##


# Together
def _read_metric(path_results, pattern=None, metrics=None):
    df = get_metrics(path_results, pattern)
    df = df.rename(columns={'metric_x':'metric', 'metric_y':'distance_metric'})[['metric', 'metric_value', 'sample']+metrics]
    return df

def rank_combos(df, metric_list, n_top=3):

    L = []
    for x in metric_list:
        L += df.loc[df['metric'] == x].groupby(metrics)['metric_value'].max().sort_values().index[-n_top:].to_list()
    
    s = pd.Series(L).value_counts()
    s = s.map(lambda x: f'{x}/{len(metric_list)}').to_frame('top 3 rank (within n metrics tested)')

    return s

##

df = pd.concat([
    _read_metric(path_results, 'dataset_df', metrics), 
    _read_metric(path_results, 'distance_metrics', metrics), 
    _read_metric(path_results, 'tree_metrics', metrics)
])

## Rank and score
good_mutational_properties = ['transitions_vs_transversions_ratio']
less_sparseness = ['median_n_assigned_char_per_clone', 'density', 'median_n_vars_per_cell']
robust_connectivity = ['corr', 'median_support_upmost_nodes']#, 'corr_distances']
good_relationship_with_GBC = ['ARI', 'NMI']#, 'median_NN_purity', 'AUCPR']#, 'freq_lineage_biased_muts']
relevant_metrics = [*robust_connectivity,*good_relationship_with_GBC]

rank_combos(df, relevant_metrics)


##


# Chosen default
cell_filter = 'filter2'
max_AD_counts = "2"
min_AD = 1
min_mean_AD_in_confident = "1.5"

(
    df.query('sample=="MDA_PT"')
    .query('cell_filter==@cell_filter and max_AD_counts==@max_AD_counts and min_AD==@min_AD and min_mean_AD_in_confident==@min_mean_AD_in_confident')
    # .query('metric in @relevant_metrics')
    # .groupby('metric')['metric_value'].describe()
)


##


# Viz
# import pickle
# with open(os.path.join(path_results, 'MDA_lung', 'job29', 'annotated_tree.pickle'), 'rb') as f:
#     tree = pickle.load(f)
# 
# _, _, muts = cut_and_annotate_tree(tree)
# 
# fig, axs = plt.subplots(1,2,figsize=(10,5))
# plot_tree(tree, ax=axs[0], orient='down', features=muts, layer='raw', feature_label_size=2)
# plot_tree(tree, ax=axs[1], orient='down', features=muts, layer='transformed', feature_label_size=2)
# fig.tight_layout()
# plt.show()
# 
# 
# fig, ax = plt.subplots(figsize=(5,5))
# # df = pd.crosstab(tree.cell_meta['GBC'], tree.cell_meta['MT_clone'], normalize=1)
# # row_sort = df.index[leaves_list(linkage(df.values))]
# # col_sort = df.columns[leaves_list(linkage(df.values.T))]
# # plot_heatmap(df.loc[row_sort,col_sort].T, ax=ax)
# df = pd.read_csv(os.path.join(path_results, 'MDA_lung', 'job29', 'evo_coupling.csv'), index_col=0)
# ax.imshow(df)
# fig.tight_layout()
# plt.show()


##