"""
Benchmarking experiment analysis
"""

import os
from itertools import chain
from mito_utils.utils import *
from mito_utils.clustering import *
from mito_utils.preprocessing import *
from mito_utils.plotting_base import *
from mito_utils.phylo_plots import *
matplotlib.use('macOSX')


##


# Get metrics
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro'
path_data = os.path.join(path_main, 'data', 'AML', 'AFMs')
path_results = os.path.join(path_main, 'results', 'AML', 'tuning')


##


# Set annot
groupings = ['pp_method', 'bin_method', 'af_confident_detection', 'min_AD', 'min_n_positive']
metric_annot = {
    'Mutation Quality' : ['n_dbSNP', 'n_REDIdb', 'transitions_vs_transversions_ratio'],
    'Association with GBC' : ['freq_lineage_biased_muts'],                               
    'Noise robustness' : ['corr'],
    'Connectivity' : ['density', 'transitivity', 'average_path_length', 'average_degree', 'proportion_largest_component'],
    'Variation' : ['n_aggregated_ct_groups', 'median_n_vars_per_cell'],                                                           
    'Yield' : ['n_cells', 'n_vars']                                                                
}  
relevant_metrics = list(chain.from_iterable([ metric_annot[k] for k in metric_annot ]))
relevant_metrics = [ f'{x}_rescaled' for x in relevant_metrics ]
weights = {
    'Mutation Quality': .1,
    'Association with GBC': .1,
    'Noise robustness' : .3,
    'Connectivity' : .0,
    'Variation' : .3,
    'Yield' : .1
}


##


# Extract
df, metrics, options = format_results(path_results)

# Score and rank, single task
n = 10
df_ranked = rank_items(df, groupings, metrics, weights, metric_annot)
df_final = pd.concat([df_ranked.head(n), df_ranked.tail(n)])
metric_type_scores = df_final.columns[df_final.columns.str.contains('score')].to_list()
df_final[groupings+metric_type_scores+relevant_metrics]

# Options of interests
sample = 'AML2'
n_top = 10
df_ranked = rank_items(df.query('sample==@sample'), 'job_id', metrics, weights, metric_annot)
df_final = df_ranked.head(n)
df_final = df_final.merge(df[['job_id']+options.to_list()], on='job_id')
metric_type_scores = df_final.columns[df_final.columns.str.contains('score')].to_list()
df_final = df_final[['job_id']+groupings+metric_type_scores+['n_cells', 'n_vars', 'transitions_vs_transversions_ratio', 'freq_lineage_biased_muts', 'corr']]




for job_id in df_final['job_id'].unique():

    afm = sc.read(os.path.join(path_data, sample, 'afm.h5ad'))
    with open(os.path.join(path_results, sample, f'{job_id}_stats.pickle'), 'rb') as f:
        d = pickle.load(f)

    afm = filter_cells(afm, cell_filter='filter2')
    afm.obs['lymph_vs_myelo'] = np.select(
        [afm.obs['aggregated_ct'].isin(['T_CD8', 'T_CD4']), afm.obs['aggregated_ct'].isin(['Early_myeloid', 'HSC_MPP', 'B_early'])],
        ['lympho', 'progenitors'], default='Mature myelo'
    ) 

    afm = filter_afm(afm, filtering=None, cells=d['cells'].to_list(), variants=d['vars'].to_list(), bin_method='vanilla', binarization_kwargs={'t_vanilla':0, 'min_AD':d['genotyping']['min_AD']})
    tree = build_tree(afm, metric='jaccard', bin_method='vanilla', solver='NJ')
    _, mut_nodes, mut_order = cut_and_annotate_tree(tree)

    fig, ax = plt.subplots(figsize=(4.8,5))
    cmaps = {
        'tumor_tme_class_new' : {'malignant':'red', 'tme':'b', 'undefined':'grey'},
        'lymph_vs_myelo' : { 'lympho' : 'green', 'progenitors' : 'orange', 'Mature myelo' : 'darkred'},
        'MT_clone' : create_palette(tree.cell_meta, 'MT_clone', sc.pl.palettes.godsnot_102)
    }
    plot_tree(tree, features=['tumor_tme_class_new', 'lymph_vs_myelo', 'MT_clone'], categorical_cmaps=cmaps, ax=ax, colorstrip_width=5)
    ax.set(title=f'n_cells: {afm.shape[0]}, n_vars {afm.shape[1]}')
    fig.tight_layout()
    fig.savefig(os.path.join(path_results, sample, f'{job_id}_tree.png'))

    fig, ax = plt.subplots(figsize=(4.8,5))
    order=leaves_list(linkage(afm.obsp['distances'].A))
    ax.imshow(afm.obsp['distances'].A[np.ix_(order, order)], cmap='inferno_r')
    fig.tight_layout()
    fig.savefig(os.path.join(path_results, sample, f'{job_id}_distances.png'), dpi=500)

    fig, ax = plt.subplots(figsize=(4,4))
    ax.plot(
        np.sum(afm[afm.obs['tumor_tme_class_new']=='tme'].layers['bin'].A==1, axis=0) / np.sum(afm.obs['tumor_tme_class_new']=='tme'),
        np.sum(afm[afm.obs['tumor_tme_class_new']=='malignant'].layers['bin'].A==1, axis=0) / np.sum(afm.obs['tumor_tme_class_new']=='malignant'),
        'ko'
    )
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set(xlim=(-0.1,1.1), ylim=(-0.1,1.1))
    format_ax(ax=ax, xlabel='TME', ylabel='Tumor')
    fig.tight_layout()
    fig.savefig(os.path.join(path_results, sample, f'{job_id}_tme_tumor.png'))

    plt.close()






##