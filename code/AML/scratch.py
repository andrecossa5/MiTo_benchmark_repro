"""
MT-SNVs tree and mutation scratch.
"""

import os
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.dimred import *
from mito_utils.kNN import *
from mito_utils.distances import *
from mito_utils.metrics import *
from mito_utils.plotting_base import *
from mito_utils.embeddings_plots import *
from mito_utils.diagnostic_plots import *
from mito_utils.phylo import *
from mito_utils.phylo_plots import *


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro'
path_data = os.path.join(path_main, 'data', 'AML')
path_results = os.path.join(path_main, 'results',  'AML')
make_folder(path_results, 'trees', overwrite=False)

# Params
sample = 'sAML1'
filtering = 'MI_TO'

# Make an annotated AF matrix for the sample
path_afm = os.path.join(path_data, 'AFMs', f'{sample}.h5ad')
path_meta = os.path.join(path_data, 'cells_meta.csv')
afm = make_AFM(path_afm, path_meta, sample=sample, cell_filter='filter2', 
               nmads=5, mean_cov_all=25, median_cov_target=20, min_perc_covered_sites=0.75, is_covered_treshold=0)
a, report = filter_AFM(
    afm, filtering=filtering, 
    filtering_kwargs={'min_median_af':.01, 'af_confident_detection': .02, 'min_n_confidently_detected':2, 'min_n_positive':5, 'min_mean_AD_in_confident':1.5 },
    max_AD_counts=2, min_cell_number=10,
    lineage_column='GBC', compute_enrichment=True, spatial_metrics=False,
    binarization_kwargs={'min_AD':2, 't_vanilla':0},
    path_dbSNP=os.path.join(path_data, 'miscellanea', 'dbSNP_MT.txt'),
    path_REDIdb=os.path.join(path_data, 'miscellanea', 'REDIdb_MT.txt'),
)

s = (a.var.loc[:,a.var.columns.str.startswith('FDR')]<=.1).sum(axis=0)
np.sum(s>0) / s.size

AD, DP, _ = get_AD_DP(a)
AD = AD.A.T
DP = DP.A.T

# (a.var.loc[:,a.var.columns.str.startswith('FDR')]<=.1).sum(axis=0)

ref_df = pd.read_csv('/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/data/MI_TO_bench/miscellanea/weng2024_mut_spectrum_ref.csv', index_col=0)
mut_profile(a, ref_df=ref_df)
plt.show()

AD, DP, _ = get_AD_DP(a)
AD = AD.A.T
DP = DP.A.T


##


tree = build_tree(AD=AD, DP=DP, meta=a.obs, var_names=a.var_names, solver='NJ',
                        metric='jaccard', bin_method='vanilla', 
                        binarization_kwargs={'t_vanilla':0.01, 't_prob':.5, 'min_AD':2})

tree, mut_nodes, mutation_order =  cut_and_annotate_tree(tree, n_clones=5)


cmaps = {
    'tumor_tme_class_new' : {'malignant':'orange', 'tme':'green'},
    'MT_clone' : create_palette(tree.cell_meta, 'MT_clone', sc.pl.palettes.vega_20)
}


fig, ax = plt.subplots(figsize=(10,10))
plot_tree(tree, ax=ax, 
          orient='down',
          colorstrip_spacing=.1, colorstrip_width=5, 
          features=['MT_clone']+mutation_order, 
          layer='transformed')
          # categorical_cmaps=cmaps,
          # internal_node_subset=mut_nodes,
          # feature_label_size=5,
          # internal_node_kwargs={'markersize':5, 'c':'darkred'})#, show_internal=True)
fig.tight_layout()
plt.show()


mut_profile(a, ref_df)
plt.show()






fig, axs = plt.subplots(1,2,figsize=(10,4.5))

plot_tree(tree, ax=axs[0], 
    colorstrip_spacing=.5, colorstrip_width=3,
    orient='down',
    features=['tumor_tme_class_new', 'MT_clone']+mutation_order, layer='raw',
    categorical_cmaps=cmaps,
    feature_label_size=8, feature_label_offset=2,
)

plot_tree(tree, ax=axs[1], 
    colorstrip_spacing=.5, colorstrip_width=3,
    orient='down',
    features=['tumor_tme_class_new', 'MT_clone']+mutation_order, layer='transformed',
    feature_label_size=8, feature_label_offset=2,
    categorical_cmaps=cmaps,
    internal_node_subset=mut_nodes,
    internal_node_kwargs={'markersize':5, 'c':'darkred'}, show_internal=True
)

fig.tight_layout()
plt.show()


##


biases_df = compute_lineage_biases(a, lineage_column='tumor_tme_class_new', target_lineage='malignant', bin_method='vanilla',
                       binarization_kwargs={'min_AD':2, 't_vanilla':.01})

biases_df.query('FDR<.1')
