"""
Distance and tree function final refactoring.
"""

import os
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.dimred import *
from mito_utils.kNN import *
from mito_utils.distances import *
from mito_utils.metrics import *
from mito_utils.plotting_base import *
from mito_utils.clustering import *
from mito_utils.embeddings_plots import *
from mito_utils.diagnostic_plots import *
from mito_utils.phylo import *
from mito_utils.phylo_plots import *


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro'
path_data = os.path.join(path_main, 'data', 'MI_TO_bench')
path_results = os.path.join(path_main, 'results',  'MI_TO_bench')
make_folder(path_results, 'discretization', overwrite=False)

# Params
sample = 'MDA_clones'
filtering = 'MI_TO'

# Make an annotated AF matrix for the sample
path_afm = os.path.join(path_data, 'AFMs', f'{sample}.h5ad')
path_meta = os.path.join(path_data, 'cells_meta.csv')
afm = make_AFM(path_afm, path_meta=None, sample=sample, cell_filter='filter1', 
               nmads=5, mean_cov_all=25, median_cov_target=20, min_perc_covered_sites=0.75, is_covered_treshold=0)
a, report = filter_AFM(
    afm, filtering=filtering, 
    filtering_kwargs={'min_median_af':.01, 'af_confident_detection': .02, 'min_n_confidently_detected':2, 'min_n_positive':5, 'min_mean_AD_in_confident':1.5 },
    max_AD_counts=2, min_cell_number=10,
    # lineage_column='GBC', compute_enrichment=True, spatial_metrics=True,
    binarization_kwargs={'min_AD':2, 't_vanilla':0},
    path_dbSNP=os.path.join(path_data, 'miscellanea', 'dbSNP_MT.txt'),
    path_REDIdb=os.path.join(path_data, 'miscellanea', 'REDIdb_MT.txt'),
)

np.median(afm.layers['coverage'].sum(axis=1))

s = (a.var.loc[:,a.var.columns.str.startswith('FDR')]<=.1).sum(axis=0)
np.sum(s>0) / s.size

AD, DP, _ = get_AD_DP(a)
AD = AD.A.T
DP = DP.A.T


##


report


# AUC

# MI_TO
# labels = a.obs['GBC']
# for t in np.linspace(.1,.9,5):
#     D = compute_distances(a=a, metric='custom_MI_TO_jaccard', bin_method='MI_TO', binarization_kwargs={'t_prob':t})
#     distance_AUPRC(D, labels)
# 
# # Vanilla
# D = compute_distances(a=a, metric='jaccard', bin_method='vanilla', binarization_kwargs={'t_vanilla':.05})
# distance_AUPRC(D, labels)
# 
# # kNN graph (fixed k, varying metric) has the
# results = {}
# k = 10
# labels = a.obs['GBC'].astype('str')
# metrics = ['custom_MI_TO_jaccard', 'jaccard']
# bin_method = ['MI_TO', 'vanilla']
# 
# for metric,bin_method in zip(metrics, bin_method):
#     D = compute_distances(a=a, metric=metric, bin_method=bin_method)
#     idx = kNN_graph(D=D, k=k, from_distances=True)[0]
#     _, _, acc_rate = kbet(idx, labels, only_score=False)
#     median_entropy = NN_entropy(idx, labels)
#     median_purity = NN_purity(idx, labels)
#     results[metric] = {
#         'kBET_rejection_rate' : 1-acc_rate, 
#         'median_NN_entropy': median_entropy, 
#         'median_NN_purity' : median_purity
#     }
# 
# # 3. What metric is more robust to noise??
# results = {}
# n_samples = 100
# n_vars_sampling = round((a.shape[1] / 100) * 80)
# metrics = ['custom_MI_TO_jaccard', 'jaccard']
# bin_method = ['MI_TO', 'vanilla']
# 
# for metric,bin_method in zip(metrics, bin_method):
#     L = []
#     for _ in range(n_samples):
#         idx = np.random.choice(range(AD.shape[1]), size=n_vars_sampling, replace=False)
#         D = compute_distances(AD=AD[:,idx], DP=DP[:,idx], metric=metric, bin_method=bin_method)
#         L.append(D.flatten())
#     results[metric] = np.mean(np.corrcoef(np.array(L)))


##


# TREES = []
# for _ in range(20):
#     tree = bootstrap_iteration(AD=AD, DP=DP, meta=a.obs, solver='NJ', metric='custom_MI_TO_jaccard', 
#                         boot_strategy='jacknife', bin_method='MI_TO', binarization_kwargs={'t_prob':.5})
#     TREES.append(tree)
# 
# tree = build_tree(AD=AD, DP=DP, meta=a.obs, solver='NJ', var_names=a.var_names,
#                   metric='custom_MI_TO_jaccard', bin_method='MI_TO', binarization_kwargs={'t_prob':.5})
# 
# tree = calculate_supports(tree, TREES)


##


# ref_df = pd.read_csv('/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/data/MI_TO_bench/miscellanea/weng2024_mut_spectrum_ref.csv', index_col=0)
# 
# for t in np.arange(1,6):
#     t = 6
#     a, report = filter_cells_and_vars(
#         afm, filtering=filtering, 
#         max_AD_counts=t, af_confident_detection=.001, min_cell_number=10,
#         lineage_column='GBC', compute_enrichment=True,
#         path_dbSNP=os.path.join(path_data, 'miscellanea', 'dbSNP_MT.txt'),
#         path_REDIdb=os.path.join(path_data, 'miscellanea', 'REDIdb_MT.txt'),
#     )
#     fig = mut_profile(a, ref_df, figsize=(6,3))
#     fig.suptitle(f'Max n ALT alleles: {t} (n MT-SNV={a.var_names.size})')
#     fig.savefig(os.path.join(path_results, f'mut_spectrum_{t}.png'), dpi=300)


##


G_smooth = call_genotypes(AD=AD, DP=DP, bin_method='MI_TO_smooth', t_prob=.95, k=5, gamma=.3, n_samples=10)
G_non_smooth = call_genotypes(AD=AD, DP=DP, bin_method='MI_TO', t_prob=.95)
G_vanilla = call_genotypes(AD=AD, DP=DP, bin_method='vanilla', t_vanilla=0.01)


# (G_smooth==1).sum(axis=0).mean()
# (G_non_smooth==1).sum(axis=0).mean()
# (G_vanilla==1).sum(axis=0).mean()
# 
# 
# idx = 3
# df = genotype_mix(AD[:,idx], DP[:,idx], t_prob=.9, t_vanilla=.05, debug=True)
# 
# df['AF'] = df['ad'] / (df['dp']+.00001)
# 
# (df['geno_vanilla']==1).sum()
# (df['geno_prob']==1).sum()
# df.loc[df['geno_vanilla'] != df['geno_prob']]

# tree_vanilla = tree_MI_TO = build_tree(AD=AD, DP=DP, meta=a.obs, var_names=a.var_names, metric='jaccard', solver='UPMGA',
#                         bin_method='vanilla', binarization_kwargs={'t_vanilla':.05})
# 
# tree_MI_TO = build_tree(AD=AD, DP=DP, meta=a.obs, var_names=a.var_names, metric='custom_MI_TO_jaccard', solver='UPMGA',
#                         bin_method='MI_TO', binarization_kwargs={'t_prob':.5})
# 
# tree_MI_TO_smooth = build_tree(AD=AD, DP=DP, meta=a.obs, var_names=a.var_names, metric='custom_MI_TO_jaccard',  solver='UPMGA',
#                        bin_method='MI_TO_smooth', binarization_kwargs={'t_prob':.9, 'k':5, 'gamma':.1, 'n_samples':10})

# fig, ax = plt.subplots(figsize=(4.5,4.5))
# 
# plot_tree(tree_MI_TO, ax=ax, 
#           colorstrip_spacing=.5,
#           orient='down',
#           features=['GBC']+tree_MI_TO.layers['transformed'].columns.to_list(), layer='transformed',
#           feature_label_size=5, feature_label_offset=2,
#           # continuous_cmaps=continuous_cmaps,
#           # feature_internal_nodes='support', internal_node_kwargs={'markersize':5}
# )
# 
# fig.tight_layout()
# plt.show()


##

report

mut_profile(a, ref_df=ref_df)
plt.show()

tree = build_tree(AD=AD, DP=DP, meta=a.obs, var_names=a.var_names, metric='jaccard', solver='UPMGA',
                      bin_method='vanilla', binarization_kwargs={'t_prob':.75, 't_vanilla':.0, 'min_AD':2})

tree, mut_nodes, mutation_order =  cut_and_annotate_tree(tree)


cmaps = {
    'GBC' : create_palette(tree.cell_meta, 'GBC', sc.pl.palettes.godsnot_102),
    'MT_clone' : create_palette(tree.cell_meta, 'MT_clone', sc.pl.palettes.default_102)
}


fig, ax = plt.subplots(figsize=(5,5))
plot_tree(tree, ax=ax)
plt.show()

fig, axs = plt.subplots(1,2,figsize=(15,8))

plot_tree(tree, ax=axs[0], 
    colorstrip_spacing=.01, colorstrip_width=2,
    orient='down',
    features=['GBC', 'MT_clone']+mutation_order, layer='raw',
    # categorical_cmaps=cmaps,
    feature_label_size=10, feature_label_offset=2,
)

plot_tree(tree, ax=axs[1],
    colorstrip_spacing=.01, colorstrip_width=2,
    orient='down',
    features=['GBC', 'MT_clone']+mutation_order, layer='transformed',
    feature_label_size=10, feature_label_offset=2,
    # categorical_cmaps=cmaps,
    internal_node_subset=mut_nodes,
    internal_node_kwargs={'markersize':5, 'c':'darkred'}, show_internal=True
)

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'PT_trees.png'))


plt.show()


from sklearn.metrics import normalized_mutual_info_score

df = tree.cell_meta.dropna()

custom_ARI(df['GBC'], df['MT_clone'])
normalized_mutual_info_score(df['GBC'], df['MT_clone'])

mut_profile(a, ref_df)
plt.show()

report


a.obs['GBC'].unique().size
