"""
Benchmark of MT-SNVs spaces, for AML samples.
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
path_data = os.path.join(path_main, 'data', 'MI_TO_bench', 'AFMs')
path_results = os.path.join(path_main, 'results', 'AML', 'tuning_maegatk')


##


# Set annot
groupings = ['pp_method', 'bin_method', 'af_confident_detection', 'min_n_confidently_detected', 'min_AD']
metric_annot = {
    'Mutation Quality' : ['n_dbSNP', 'n_REDIdb', 'transitions_vs_transversions_ratio'],
    # 'Association with GBC' : ['freq_lineage_biased_muts', 'AUPRC', 'ARI', 'NMI'],                               
    'Tree structure' : ['corr', 'mean_CI'],
    'Connectedness' : ['density', 'transitivity', 'average_path_length', 'average_degree', 'proportion_largest_component'],
    'Variation' : ['genomes_redundancy', 'median_n_vars_per_cell'],                                                           
    'Yield' : ['n_GBC_groups', 'n_cells', 'n_vars']                                                                
}  
relevant_metrics = list(chain.from_iterable([ metric_annot[k] for k in metric_annot ]))
relevant_metrics = [ f'{x}_rescaled' for x in relevant_metrics ]
weights = {
    'Mutation Quality': .2,
    # 'Association with GBC': .4,
    'Tree structure' : .5,
    'Connectedness' : .0,
    'Variation' : 0,
    'Yield' : .3
}


##


# Extract
# df, metrics, options = format_results(path_results)
df = pd.read_csv(os.path.join(path_results, 'df_tuning.csv'), index_col=0)
metrics = pd.read_csv(os.path.join(path_results, 'metrics.csv')).iloc[:,0].to_list()
metrics = [ x for x in metrics if x != 'median_target/untarget_coverage_logratio']
options = pd.read_csv(os.path.join(path_results, 'options.csv')).iloc[:,0].to_list()
options += ['pp_method']
df = df.drop(columns=['median_target/untarget_coverage_logratio'])

# One sample/task
sample = 'sAML1'
df = df.query('sample==@sample')

# Explore ...
# ...

# Score and rank, single task
# n = 5
# df_ranked = rank_items(df, groupings, metrics, weights, metric_annot)
# df_final = pd.concat([df_ranked.head(n), df_ranked.tail(n)])
# metric_type_scores = df_final.columns[df_final.columns.str.contains('score')].to_list()
# df_final[groupings+['ARI', 'NMI', 'AUPRC', 'corr', 'n_cells', 'n_vars', 'n_GBC_groups']].head(5).describe()
# df_final[metric_type_scores]
# df_final[groupings+metric_type_scores+relevant_metrics].to_csv(os.path.join(path_results, 'grouped.csv'))

# Options of interests
df_selected = (
    df.query('corr>.5 and n_cells>1000 and n_vars>10') 
    [['job_id', 'pp_method', 'bin_method', 'af_confident_detection', 'min_AD', 'corr', 'freq_lineage_biased_muts', 'n_cells', 'n_vars', 'mean_CI']]
)
df_selected
df_selected['cat'] = pd.cut(df_selected['n_vars'], bins=5)
df_selected.groupby('cat').size() 
df_selected = df_selected.groupby('cat').apply(lambda x: x.sort_values('corr', ascending=False).head(2)).reset_index(drop=True)
df_selected

# Write out
# L = []
# for i in range(df_selected.shape[0]):
#     l = [sample, os.path.join(path_data, df_selected['pp_method'].values[i], sample, 'afm.h5ad'), df_selected['job_id'].values[i], "None"]
#     L.append(l)
# pd.DataFrame(L, columns=['sample', 'ch_matrix', 'job_id', 'cell_file']).set_index('sample').to_csv(os.path.join(path_results, 'final_jobs.csv'))


##



# FILTER USED!!, THESIS
# MDA_clones: df.query('AUPRC>.5 and corr>.6 and n_cells>300 and n_GBC_groups==7 and n_vars>10')    
# MDA_PT: df.query('AUPRC>.3 and corr>.5 and n_cells>1000 and n_GBC_groups>30 and n_vars>10')       # q=10, n=1
# MDA_lung: df.query('AUPRC>.5 and corr>.5 and n_cells>1000 and n_GBC_groups>8 and n_vars>10')      # q=10, n=1

# ORIGINAL
# weights = {
#     'Mutation Quality': .1,
#     'Association with GBC': .4,
#     'Tree structure' : .2,
#     'Connectedness' : .0,
#     'Variation' : .0,
#     'Yield' : .3
# }
# NOW
# weights = {
#     'Mutation Quality': .1,
#     'Association with GBC': .5,
#     'Tree structure' : .1,
#     'Connectedness' : .0,
#     'Variation' : 0,
#     'Yield' : .3
# }


##