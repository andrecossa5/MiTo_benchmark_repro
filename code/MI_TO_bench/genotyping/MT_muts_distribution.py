"""
MT-SNVs AF distribution and fit to known distribution families. 
"""

import os
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.genotyping import *
from mito_utils.dimred import *
from mito_utils.plotting_base import *
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
afm = read_one_sample(path_afm, path_meta, sample=sample, nmads=5, mean_coverage=25)

# Read and filter
_, a = filter_cells_and_vars(afm, filtering=filtering, max_AD_counts=3, af_confident_detection=0, min_cell_number=10)
AD, DP, _ = get_AD_DP(a)
AD = AD.A.T
DP = DP.A.T

# VAF spectrum
fig, ax = plt.subplots(figsize=(4.5,4.5))
vars_AF_dist(a, ax=ax, color='darkred', linewidth=1, linestyle='-')
fig.tight_layout()
plt.show()

(a.X>0).sum(axis=1)                 # n vars (per cell)
(a.X>0).sum(axis=0) / a.shape[0]    # Prop +cells (per var)
a.X.mean(axis=0)                    # Mean AF vars
AD.mean(axis=0)                     # Mean AD vars
DP.mean(axis=0)                     # Mean DP vars

fig, ax = plt.subplots(figsize=(4.5,4.5))
ax.plot((a.X>0).sum(axis=0), np.ma.masked_less_equal(AD, 0).mean(axis=0), 'ko')
format_ax(ax=ax, xlabel='n +cells', ylabel='Mean n of UMIs for ALT in +cells')
fig.tight_layout()
plt.show()


##


# Benchmark
geno = genotype_MI_TO(a, debug=True)
geno.sum(axis=1)

for t in np.linspace(.5,.9,5):
    labels = a.obs['GBC']
    evaluate_metric_with_gt(
        a, metric='custom_MI_TO_jaccard', labels=labels, bin_method='MI_TO', 
        discretization_kwargs={'t_prob':t}
    )
    evaluate_metric_with_gt(
        a, metric='jaccard', labels=labels, bin_method='vanilla'
    )


##


X, dimnames = reduce_dimensions(a, method='UMAP', metric='cosine', n_comps=2)
df_ = pd.DataFrame(X, columns=dimnames, index=a.obs_names).join(a.obs)

fig, ax = plt.subplots(figsize=(4.5,4.5))
draw_embeddings(df_, cat='GBC', ax=ax, legend_kwargs={'loc':'upper right', 'bbox_to_anchor':(1,1)})
fig.tight_layout()
plt.show()


tree = build_tree(a, solver='NJ')
tree.cell_meta['GBC'] = a.obs['GBC'].astype(str)
tree.cell_meta = tree.cell_meta.join(pd.DataFrame(a.X, index=a.obs_names, columns=a.var_names))

variants = get_supporting_muts(tree, a)

fig, ax = plt.subplots(figsize=(7,7))
plot_tree(tree, ax=ax, meta=['GBC'],#+variants, 
          colorstrip_width=2, colorstrip_spacing=.2, orient=90, vmin_annot=0, vmax_annot=.05)
fig.tight_layout()
plt.show()

tree.n_character

help(tree.remove_leaves_and_prune_lineages)




cs.tl.compute_expansion_pvalues(tree, min_clade_size=(0.5 * tree.n_cell), min_depth=1)
# this specifies a p-value for identifying expansions unlikely to have occurred
# in a neutral model of evolution
probability_threshold = 0.1

expanding_nodes = []
for node in tree.depth_first_traverse_nodes():
    if tree.get_attribute(node, "expansion_pvalue") < probability_threshold:
        expanding_nodes.append(node)






tree.collapse_mutationless_edges(True)

dir(tree)

len(tree.internal_nodes)

CI(tree)
RI(tree)


[ len(v) for k,v in get_clades(tree).items() ]

cells_ = list(list(get_clades(tree).values())[1])
list(get_clades(tree).keys())[1]

tree.subset_clade('cassiopeia_internal_nodeba35136b13cf23c5f60364a2')


tree



fig, ax = plt.subplots(figsize=(7,7))
plot_tree(tree, ax=ax, meta=['GBC'],#+variants, 
          colorstrip_width=2, colorstrip_spacing=.2, orient=90, vmin_annot=0, vmax_annot=.05)
fig.tight_layout()
plt.show()




cells = t_removed.cell_meta.query('GBC=="TGCAGTTTTGGTGCTCTA"').index



t_ = t_removed.copy()
T = t_.get_tree_topology()

for cell in cells:
    T.remove_node(cell)

cells = t_.character_matrix.index[t_.character_matrix.index.isin(T.nodes)]
t_removed_2 = CassiopeiaTree(
    tree=T, 
    character_matrix=t_.character_matrix.loc[cells],
    cell_meta=t_.cell_meta.loc[cells,:]
)

t_removed_2.leaves
t_removed_2.cell_meta


fig, ax = plt.subplots(figsize=(7,7))
plot_tree(t_removed, ax=ax, meta=['GBC'],#+variants, 
          colorstrip_width=2, colorstrip_spacing=.2, orient=90, vmin_annot=0, vmax_annot=.05)
fig.tight_layout()
plt.show()

t_removed.cell_meta['GBC'].value_counts()