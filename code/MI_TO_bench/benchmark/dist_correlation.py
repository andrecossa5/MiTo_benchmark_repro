"""
Contribute of individual options to individual metrics.
"""

import os
from itertools import product
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.plotting_base import *
matplotlib.use('macOSX')


##


# Get metrics
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro'
path_data = os.path.join(path_main, 'results', 'MI_TO_bench', 'benchmark', 'tuning')
path_results = os.path.join(path_main, 'results', 'MI_TO_bench', 'benchmark')

# Main metrics and options df
# df, metrics, options = format_results(path_data)
df = pd.read_csv(os.path.join(path_data, 'main_df.csv'), index_col=0)
metrics = pd.read_csv(os.path.join(path_data, 'metrics.csv')).iloc[:,0].to_list()
metrics = [ x for x in metrics if x != 'median_target/untarget_coverage_logratio']
options = pd.read_csv(os.path.join(path_data, 'options.csv')).iloc[:,0].to_list()
df = df.drop(columns=['median_target/untarget_coverage_logratio'])


##


# Overall
weights = {
    'Mutation Quality': .1,
    'Association with GBC': .4,
    'Tree structure' : .2,
    'Connectedness' : .0,
    'Variation' : .0,
    'Yield' : .3
}
metric_annot = {
    'Mutation Quality' : ['n_dbSNP', 'n_REDIdb', 'transitions_vs_transversions_ratio'],
    'Association with GBC' : ['freq_lineage_biased_muts',  'AUPRC', 'ARI', 'NMI'],                               
    'Tree structure' : ['corr'],
    'Connectedness' : ['density', 'transitivity', 'average_path_length', 'average_degree', 'proportion_largest_component'],
    'Variation' : ['genomes_redundancy', 'median_n_vars_per_cell', 'n_vars'],                                                           
    'Yield' : ['n_GBC_groups', 'n_cells']                                                                
}  
groupings = ['pp_method', 'bin_method', 'af_confident_detection', 'min_AD', 'min_n_positive']

df_ = rank_items(df, groupings, metrics, weights, metric_annot)
df_.columns


##


fig, axs = plt.subplots(1,4,figsize=(13,3.5))

metrics = ['freq_lineage_biased_muts', 'AUPRC', 'ARI', 'NMI']

for i,metric in enumerate(metrics):
    ax = axs.ravel()[i]
    sns.regplot(df_['corr'], df_[metric], ax=ax, scatter=False)
    ax.plot(df_['corr'], df_[metric], color='lightgrey', 
            marker='o', linestyle='', markersize=5, markeredgecolor='k')
    corr = np.corrcoef(df_['corr'].values, df_[metric].values)[0,1]
    ax.set(title=f'Pearson\'s r: {corr:.2f}')

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'dist_corr_and_other_GBC.png'), dpi=500)


##