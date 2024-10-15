"""
Benchmarking experiment analysis
"""

import os
from itertools import chain
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.plotting_base import *
matplotlib.use('macOSX')


##


# Get metrics
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro'
path_data = os.path.join(path_main, 'results', 'MI_TO_bench', 'phylo_inference')
path_results = os.path.join(path_main, 'results', 'MI_TO_bench', 'benchmark')


##


# Set annot
groupings = ['pp_method', 'bin_method', 'af_confident_detection', 'min_AD', 'min_n_positive']
metric_annot = {
    'Mutation Quality' : ['median_target/untarget_coverage_logratio', 'n_dbSNP', 'n_REDIdb', 'transitions_vs_transversions_ratio'],
    'Association with GBC' : ['freq_lineage_biased_muts',  'AUPRC', 'ARI', 'NMI'],                               
    'Noise robustness' : ['corr'],
    'Connectivity' : ['density', 'transitivity', 'average_path_length', 'average_degree', 'proportion_largest_component'],
    'Variation' : ['genomes_redundancy', 'median_n_vars_per_cell'],                                                           
    'Yield' : ['n_GBC_groups', 'n_cells']                                                                
}  
relevant_metrics = list(chain.from_iterable([ metric_annot[k] for k in metric_annot ]))
relevant_metrics = [ f'{x}_rescaled' for x in relevant_metrics ]
weights = {
    'Mutation Quality': 0.1,
    'Association with GBC': 0.5,
    'Noise robustness' : .3,
    'Connectivity' : .0,
    'Variation' : .0,
    'Yield' : .2
}


##


# Extract
df, metrics, options = format_results(path_data)
df.to_csv(os.path.join(path_results, 'main_benchmarking_df.csv'))

# One sample/task
df = df.query('sample=="MDA_lung"')

# Score and rank, single task
n = 5
df_ranked = rank_items(df, groupings, metrics, weights, metric_annot)
df_final = pd.concat([df_ranked.head(n), df_ranked.tail(n)])
metric_type_scores = df_final.columns[df_final.columns.str.contains('score')].to_list()
df_final[groupings+metric_type_scores+relevant_metrics].to_csv(os.path.join(path_results, 'grouped.csv'))
# df_final[groupings+metric_type_scores+relevant_metrics].columns

# Options of interests
n_top = 5
df_ranked = rank_items(df, 'job_id', metrics, weights, metric_annot)
df_final = pd.concat([df_ranked.head(n), df_ranked.tail(n)])
df_final = df_final.merge(df[['job_id']+options.to_list()], on='job_id')
metric_type_scores = df_final.columns[df_final.columns.str.contains('score')].to_list()
df_final[groupings+metric_type_scores].to_csv(os.path.join(path_results, 'single_jobs.csv'))

# df_final.iloc[0,:][options]


##


# Viz option impact single metrics








