"""
Coding game.
"""

import os
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.phylo import *


sample ='MDA_PT'
path_afm = f'/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/data/MI_TO_bench/AFMs/mito_preprocessing/{sample}/afm.h5ad'
path_dbSNP = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/data/MI_TO_bench/miscellanea/dbSNP_MT.txt'
path_REDIdb = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/data/MI_TO_bench/miscellanea/REDIdb_MT.txt'


afm = sc.read(path_afm)
afm = filter_cells(afm, cell_filter='filter2')
afm, tree = filter_afm(
    afm,
    min_cell_number=10,
    lineage_column='GBC',
    filtering_kwargs={
        'min_cov' : 10,
        'min_var_quality' : 30,
        'min_frac_negative' : .2,
        'min_n_positive' : 2,
        'af_confident_detection' : .01,
        'min_n_confidently_detected' : 2,
        'min_mean_AD_in_positives' : 1.25,       # 1.25,
        'min_mean_DP_in_positives' : 25
    },
    binarization_kwargs={
        't_prob':.75, 't_vanilla':.001, 'min_AD':1, 'min_cell_prevalence':.01
    },
    bin_method='MI_TO',
    tree_kwargs={'metric':'custom_MI_TO_jaccard', 'solver':'NJ', 'ncores' : 8},
    path_dbSNP=path_dbSNP, 
    path_REDIdb=path_REDIdb,
    spatial_metrics=True,
    compute_enrichment=False,
    max_AD_counts=2,
    return_tree=True
)

# tree, _, _ = cut_and_annotate_tree(tree)
# stats = { k:v for k,v in afm.uns.items() }
# stats['n_MT_clone'] = tree.cell_meta['MT_clone'].nunique()
# stats['corr_dist'] = calculate_corr_distances(tree)



