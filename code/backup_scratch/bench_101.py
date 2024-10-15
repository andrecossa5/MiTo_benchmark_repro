"""
Scratch functions to create trees from other preprocessing pipelines.
"""

import os
from mito_utils.preprocessing import *
from mito_utils.metrics import *
from mito_utils.clustering import *
from mito_utils.phylo import *
from sklearn.metrics import normalized_mutual_info_score


# Path
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/'
path_data = os.path.join(path_main, 'data', 'MI_TO_bench', 'AFMs')
sample = 'MDA_PT'
lineage_column = 'GBC'
# /Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/data/MI_TO_bench/miscellanea/REDIdb_MT.txt
# /Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/data/MI_TO_bench/miscellanea/REDIdb_MT.txt


##


# Load AFM


# Binarize and select cells


X = AD / (DP+.0000001)
X.index = X.index.map(lambda x: f'{x}_{sample}')
cells = mytree.leaves
t = .1

X_bin = pd.DataFrame(np.where(X.loc[cells,:]>t,1,0), index=cells, columns=X.columns)
cells = ((X_bin>0).sum(axis=1)>1).loc[lambda x: x].index
X_bin = X_bin.loc[cells,:]

# Tree
D = pairwise_distances(X_bin.values, metric='jaccard')
D = pd.DataFrame(D, index=cells, columns=cells)
tree = CassiopeiaTree(character_matrix=X_bin, dissimilarity_map=D, cell_meta=mytree.cell_meta.loc[cells])
solver = cs.solver.NeighborJoiningSolver(add_root=True)
solver.solve(tree)
tree.layers['raw'] = X.loc[cells,:]
tree.layers['transformed'] = X_bin
tree, _, _ = cut_and_annotate_tree(tree)

# Scoring
stats = { k:v for k,v in afm.uns.items() }
stats['n_MT_clone'] = tree.cell_meta['MT_clone'].nunique()
stats['corr_dist'] = calculate_corr_distances(tree)

if lineage_column is not None:
    lineage_metrics = {}
    lineage_metrics[f'n_{lineage_column}_groups'] = tree.cell_meta[lineage_column].nunique()
    lineage_metrics['AUPRC'] = distance_AUPRC(afm.obsp['distances'].A, afm.obs[lineage_column])
    lineage_metrics['ARI'] = custom_ARI(tree.cell_meta[lineage_column], tree.cell_meta['MT_clone'])
    lineage_metrics['NMI'] = normalized_mutual_info_score(tree.cell_meta.dropna()[lineage_column], tree.cell_meta.dropna()['MT_clone'])
    stats['lineage_metrics'] = lineage_metrics

# Save
job_id = generate_job_id()
with open(f'tuning{job_id}_stats.pickle', 'wb') as f:
    pickle.dump(stats, f)


##












