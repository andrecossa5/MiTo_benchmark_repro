"""
Scratch functions to create trees from other preprocessing pipelines.
"""

import os
from mito_utils.preprocessing import *
from mito_utils.diagnostic_plots import *
from mito_utils.metrics import *
from mito_utils.clustering import *
from mito_utils.distances import *
from mito_utils.phylo import *
from mito_utils.phylo_plots import *
from mito_utils.plotting_base import *
from mito_utils.dimred import *
from mito_utils.embeddings_plots import *


# Path
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/'
path_data = os.path.join(path_main, 'data', 'MI_TO_bench', 'AFMs')
sample = 'MDA_PT'


##


# Mine top choice
with open(os.path.join(path_main, 'results', 'MI_TO_bench', 'phylo_inference', 'MDA_PT', 'job53', 'annotated_tree.pickle'), 'rb') as f:
    mytree = pickle.load(f)
cells = mytree.leaves

# Samtools/freebayes
at = pd.read_csv(os.path.join(path_data, 'freebayes', sample, 'allele_table.csv.gz'), index_col=0)
at['mut'] = at['POS'].astype(str) + at['REF'] + '>' + at['ALT']
AD = at.pivot(index='cell', columns='mut', values='AD').fillna(0)
DP = at.pivot(index='cell', columns='mut', values='DP').fillna(0)

# Maegatk
maegatk_d = { 
    '>'.join(x.split('.')[:-1]) : pd.read_csv(os.path.join(path_data, 'maegatk', sample, 'tables',x), index_col=0)
    for x in os.listdir(os.path.join(path_data, 'maegatk', sample, 'tables')) if x.startswith('n') 
}
maegatk_d.keys()
X = maegatk_d['n10>5'] / 100






X = AD / (DP+.0000001)
X.index = X.index.map(lambda x: f'{x}_{sample}')
cells = mytree.leaves
t = .1

X_bin = pd.DataFrame(np.where(X.loc[cells,:]>t,1,0), index=cells, columns=X.columns)
cells = ((X_bin>0).sum(axis=1)>1).loc[lambda x: x].index
X_bin = X_bin.loc[cells,:]
D = pairwise_distances(X_bin.values, metric='jaccard')
D = pd.DataFrame(D, index=cells, columns=cells)
tree = CassiopeiaTree(character_matrix=X_bin, dissimilarity_map=D, cell_meta=mytree.cell_meta.loc[cells])
solver = cs.solver.NeighborJoiningSolver(add_root=True)
solver.solve(tree)
tree.layers['raw'] = X.loc[cells,:]
tree.layers['transformed'] = X_bin
tree, _, _ = cut_and_annotate_tree(tree)

tree.collapse_mutationless_edges(True)
mytree.collapse_mutationless_edges(True)

fig, axs = plt.subplots(1,2,figsize=(10,5))
treecmaps = {
    'GBC' : create_palette(tree.cell_meta, 'GBC', sc.pl.palettes.default_102),
    'MT_clone' : create_palette(tree.cell_meta.dropna(), 'MT_clone', sc.pl.palettes.godsnot_102)
}
plot_tree(tree, ax=axs[0], features=['GBC', 'MT_clone'], colorstrip_width=1, categorical_cmaps=treecmaps)
mycmaps = {
    'GBC' : create_palette(mytree.cell_meta, 'GBC', sc.pl.palettes.default_102),
    'MT_clone' : create_palette(mytree.cell_meta, 'MT_clone', sc.pl.palettes.godsnot_102)
}
plot_tree(mytree, ax=axs[1], features=['GBC', 'MT_clone'], colorstrip_width=1, categorical_cmaps=mycmaps)
fig.tight_layout()
plt.show()



# harmonize_colors(tree.cell_meta, 'GBC', 'MT_clone')


from sklearn.metrics import *

distance_AUPRC(D.values, labels=tree.cell_meta['GBC'])
distance_AUPRC(mytree.get_dissimilarity_map().values, labels=mytree.cell_meta['GBC'])
calculate_corr_distances(tree)
calculate_corr_distances(mytree)
custom_ARI(tree.cell_meta['GBC'], tree.cell_meta['MT_clone'])
custom_ARI(mytree.cell_meta['GBC'], mytree.cell_meta['MT_clone'])
normalized_mutual_info_score(tree.cell_meta.dropna()['GBC'], tree.cell_meta.dropna()['MT_clone'])
normalized_mutual_info_score(mytree.cell_meta['GBC'], mytree.cell_meta['MT_clone'])

ref_df = pd.read_csv(os.path.join(path_main, 'data/MI_TO_bench/miscellanea/weng2024_mut_spectrum_ref.csv'), index_col=0)

tree_list = tree.character_matrix.columns
fig = mut_profile(mut_list=tree_list, ref_df=ref_df)
fig.savefig(os.path.join(path_results, 'mut_profiles', 'VG_Lareu.png'), dpi=1000)

my_list = mytree.character_matrix.columns
fig = mut_profile(mut_list=my_list, ref_df=ref_df)
fig.savefig(os.path.join(path_results, 'mut_profiles', 'mito_prep.png'), dpi=1000)
plt.show()





fig, axs = plt.subplots(1,2,figsize=(10,5))
embs = reduce_dimensions(X=X_bin, metric='cosine', n_comps=2)
df = tree.cell_meta.join(embs)
# axs[0].plot(df['UMAP1'], df['UMAP2'], 'ko', markersize=1)
draw_embeddings(df, cat='GBC', ax=axs[0])

embs = reduce_dimensions(X=mytree.layers['raw'], metric='cosine', n_comps=2)
df = mytree.cell_meta.join(embs)
# axs[1].plot(df['UMAP1'], df['UMAP2'], 'ko', markersize=1)
draw_embeddings(df, cat='GBC', ax=axs[1])
plt.show()

fig, axs = plt.subplots(1,2,figsize=(10,5))

leaves_list(linkage(pairwise_distances(D)))
axs[0].imshow(pd.crosstab(tree.cell_meta['GBC'], tree.cell_meta['MT_clone'], normalize=1))
axs[1].imshow(pd.crosstab(mytree.cell_meta['GBC'], mytree.cell_meta['MT_clone'], normalize=1))
plt.show()

mycross = pd.crosstab(mytree.cell_meta['GBC'], mytree.cell_meta['MT_clone'], normalize=True)
treecross = pd.crosstab(tree.cell_meta['GBC'], tree.cell_meta['MT_clone'], normalize=True)

-np.sum(mycross.values.flatten() * np.log10(mycross.values.flatten()+.00001))
-np.sum(treecross.values.flatten() * np.log10(treecross.values.flatten()+.00001))


















