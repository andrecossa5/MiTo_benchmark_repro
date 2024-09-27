"""
Discretization. NB and tresholds.
"""

import os
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.genotyping import *
from mito_utils.genotyping import _genotype_mix
from mito_utils.dimred import *
from mito_utils.kNN import *
from mito_utils.metrics import *
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
a, report = filter_cells_and_vars(
    afm, filtering=filtering, 
    max_AD_counts=3, af_confident_detection=.01, min_cell_number=10,
    lineage_column='GBC', compute_enrichment=True,
    path_dbSNP=os.path.join(path_data, 'miscellanea', 'dbSNP_MT.txt'),
    path_REDIdb=os.path.join(path_data, 'miscellanea', 'REDIdb_MT.txt'),
)

(a.var.loc[:,a.var.columns.str.startswith('FDR')]<=.1).sum(axis=0)

AD, DP, _ = get_AD_DP(a)
AD = AD.A.T
DP = DP.A.T



## ...


# tree = build_tree(a, solver='NJ')
# tree.cell_meta['GBC'] = a.obs['GBC'].astype(str)
# tree.cell_meta = tree.cell_meta.join(pd.DataFrame(a.X, index=a.obs_names, columns=a.var_names))

variants = get_supporting_muts(tree, a)

fig, ax = plt.subplots(figsize=(7,7))
plot_tree(tree, ax=ax, meta=['GBC'],#+variants, 
          colorstrip_width=2, colorstrip_spacing=.2, orient=90, vmin_annot=0, vmax_annot=.05)
fig.tight_layout()
plt.show()


len(tree.internal_nodes)
tree.collapse_mutationless_edges(True)
CI(tree)
RI(tree)


# cs.tl.compute_expansion_pvalues(tree, min_clade_size=(0.5 * tree.n_cell), min_depth=1)
# # this specifies a p-value for identifying expansions unlikely to have occurred
# # in a neutral model of evolution
# probability_threshold = 0.1
# 
# expanding_nodes = []
# for node in tree.depth_first_traverse_nodes():
#     if tree.get_attribute(node, "expansion_pvalue") < probability_threshold:
#         expanding_nodes.append(node)




# [ len(v) for k,v in get_clades(tree).items() ]
# 
# cells_ = list(list(get_clades(tree).values())[1])
# list(get_clades(tree).keys())[1]
# 
# tree.subset_clade('cassiopeia_internal_nodeba35136b13cf23c5f60364a2')
# 
# 
# tree



# fig, ax = plt.subplots(figsize=(7,7))
# plot_tree(tree, ax=ax, meta=['GBC'],#+variants, 
#           colorstrip_width=2, colorstrip_spacing=.2, orient=90, vmin_annot=0, vmax_annot=.05)
# fig.tight_layout()
# plt.show()
# 
# 
# cells = t_removed.cell_meta.query('GBC=="TGCAGTTTTGGTGCTCTA"').index
# 
# 
# 
# t_ = t_removed.copy()
# T = t_.get_tree_topology()
# 
# for cell in cells:
#     T.remove_node(cell)
# 
# cells = t_.character_matrix.index[t_.character_matrix.index.isin(T.nodes)]
# t_removed_2 = CassiopeiaTree(
#     tree=T, 
#     character_matrix=t_.character_matrix.loc[cells],
#     cell_meta=t_.cell_meta.loc[cells,:]
# )
# 
# t_removed_2.leaves
# t_removed_2.cell_meta
# 
# 
# fig, ax = plt.subplots(figsize=(7,7))
# plot_tree(t_removed, ax=ax, meta=['GBC'],#+variants, 
#           colorstrip_width=2, colorstrip_spacing=.2, orient=90, vmin_annot=0, vmax_annot=.05)
# fig.tight_layout()
# plt.show()
# 
# t_removed.cell_meta['GBC'].value_counts()