"""
Benchmark of MT-SNVs spaces.
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
# path_results = os.path.join(path_main, 'results', 'MI_TO_bench', 'benchmark', 'tuning', 'last_run_for_thesis')
path_results = os.path.join(path_main, 'results', 'MI_TO_bench', 'benchmark_final')


##


# Set annot
# groupings = ['pp_method', 'bin_method', 'af_confident_detection', 'min_n_confidently_detected', 'min_AD']
groupings = 'job_id'
metric_annot = {
    'Mutation Quality' : ['n_dbSNP', 'n_REDIdb', 'transitions_vs_transversions_ratio'],
    'Association with GBC' : ['freq_lineage_biased_muts', 'AUPRC', 'ARI', 'NMI'],                               
    'Tree structure' : ['corr', 'mean_CI'],
    'Connectedness' : ['density', 'transitivity', 'average_path_length', 'average_degree', 'proportion_largest_component'],
    'Variation' : ['genomes_redundancy', 'median_n_vars_per_cell'],                                                           
    'Yield' : ['n_GBC_groups', 'n_cells', 'n_vars']                                                                
}  
relevant_metrics = list(chain.from_iterable([ metric_annot[k] for k in metric_annot ]))
relevant_metrics = [ f'{x}_rescaled' for x in relevant_metrics ]
weights = {
    'Mutation Quality': .1,
    'Association with GBC': .4,
    'Tree structure' : .1,
    'Connectedness' : .0,
    'Variation' : 0,
    'Yield' : .4
}


##


# Extract
# df, metrics, options = format_results(path_data)
df = pd.read_csv(os.path.join(path_results, 'tuning_df.csv'), index_col=0)
metrics = pd.read_csv(os.path.join(path_results, 'metrics.csv')).iloc[:,0].to_list()
metrics = [ x for x in metrics if x != 'median_target/untarget_coverage_logratio']
options = pd.read_csv(os.path.join(path_results, 'options.csv')).iloc[:,0].to_list()
options += ['pp_method']
df = df.drop(columns=['median_target/untarget_coverage_logratio'])

# Explore
df.groupby('sample')[['corr', 'n_cells', 'n_GBC_groups', 'n_vars', 'ARI', 'NMI']].describe().T

# One sample/task
sample = 'MDA_lung'
df = df.query('sample==@sample')
df['perc_unassigned'] = df['unassigned'] / df['n_cells']
df['delta_GBC_groups'] = np.abs(df['n_GBC_groups'] - df['n MiTo clone'])
(
    df.query('n_cells>=1000 and n_GBC_groups>=9 and n_vars>10')
    [['ARI', 'corr', 'NMI', 'AUPRC', 'n_cells', 'n_vars', 'n_GBC_groups', 'n MiTo clone', 'perc_unassigned', 'delta_GBC_groups']]
)

# df.sort_values('ARI', ascending=False).head(50)[['ARI', 'corr', 'NMI', 'AUPRC', 'n_cells', 'n_vars']]
# df.query('AUPRC>.25 and corr>.25 and n_cells>=850 and n_GBC_groups>=9 and n_vars>10')[['ARI', 'corr', 'NMI', 'AUPRC', 'n_cells', 'n_vars']]

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
    df.query('n_cells>=1000 and n_GBC_groups>=9 and n_vars>10')
    [[
        'job_id', 'pp_method', 'bin_method', 'af_confident_detection',
        'ARI', 'corr', 'NMI', 'AUPRC', 'n_cells', 'n_vars', 'n_GBC_groups', 'n MiTo clone', 'perc_unassigned', 'delta_GBC_groups'
    ]]
    .sort_values('ARI', ascending=False)
    .head(5)
)

# df_selected['cat'] = pd.cut(df_selected['n_vars'], bins=3)
# df_selected.groupby('cat').size() 
# df_selected = df_selected.groupby('cat').apply(lambda x: x.sort_values('ARI', ascending=False).head(2)).reset_index(drop=True)
# df_selected[['NMI', 'ARI']].describe()

# Write out
path_ = '/data/cossa2_dare/MI_TO_bench/data/AFMs'
# L = []
for i in range(df_selected.shape[0]):
    l = [df_selected['job_id'].values[i], sample, os.path.join(path_, df_selected['pp_method'].values[i], sample, 'afm.h5ad')]
    L.append(l)

(
    pd.DataFrame(L, columns=['job_id', 'sample', 'ch_matrix'])
    .set_index('job_id')
    .to_csv(os.path.join(path_results, 'final_jobs.csv'))
)

##


# ...





#### REMINDERs


# FILTER USED!! MiTo Tree annotator final refactoring, AFMs
# 'AUPRC>=.25 and corr>=.25 and n_cells>=250(850) and n_GBC_groups>=6(8,30) and n_vars>=10'), ranked by ARI

# FILTER USED!!, All benchmark mito_preprocessing and maegatk. Before MiTo Tree annotator final refactoring
# AC_AC_PTs_1: df.query('corr>.5 and n_cells>1000 and n_GBC_groups==9 and n_vars>10')

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