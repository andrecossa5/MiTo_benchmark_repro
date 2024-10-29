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


# Annotate
metrics_ = [
    'density_rescaled', 'average_path_length_rescaled', 
    'average_degree_rescaled', 'proportion_largest_component_rescaled', 'transitivity_rescaled'
]
metric_names_ = {
    'density_rescaled':'Density', 'average_path_length_rescaled' : 'Mean path length',
    'average_degree_rescaled':'Mean degree', 
    'proportion_largest_component_rescaled':'%cells in large\nconnected component',
    'transitivity_rescaled' : 'Transitivity'
}


fig, axs = plt.subplots(1,5,figsize=(12,3), sharey=True)

for i,metric in enumerate(metrics_):
    ax = axs.ravel()[i]
    sns.regplot(df_[metric], df_['Association with GBC score'], ax=ax, scatter=False)
    ax.plot(df_[metric], df_['Association with GBC score'], color='lightgrey', 
            marker='o', linestyle='', markersize=5, markeredgecolor='k')
    name = metric_names_[metric] if metric in metric_names_ else metric
    ylabel = 'Association with GBC score' if i==0 else ''
    corr = np.corrcoef(df_[metric], df_['Association with GBC score'])[0,1]
    format_ax(ax=ax, xlabel=name, ylabel=ylabel, title=f'Pearson\'s r: {corr:.2f}', reduced_spines=True)

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'connectivity_vs_GBC.png'), dpi=500)


##



fig, axs = plt.subplots(3,5,figsize=(12,7.5), sharey=True)

# Annotate
columns_ = [
    'min_n_positive', 'af_confident_detection', 'min_mean_AD_in_positives', 'min_AD', 'median_n_vars_per_cell'
]
col_names_ = {
    'min_n_positive':'n +cells', 'af_confident_detection' : 'AF confident detection',
    'min_mean_AD_in_positives' : 'Mean n ALT UMIs in +cells', 
    'min_AD' : 'n ALT UMIs\nfor genotyping',
    'median_n_vars_per_cell' : 'Median n MT-SNVs per cell'
}
GBC_cols = ['NMI', 'ARI', 'AUPRC']

df_ = df[['pp_method', 'bin_method']+columns_+GBC_cols].dropna().groupby(['pp_method', 'bin_method']+columns_)[GBC_cols].mean().reset_index()

for i,(metric,col) in enumerate(list(product(GBC_cols, columns_))):
    ax = axs.ravel()[i]
    if col == 'median_n_vars_per_cell':
        sns.regplot(df_[col], df_[metric], ax=ax, scatter=False)
        ax.plot(df_[col], df_[metric], color='lightgrey', marker='o', linestyle='', markersize=5, markeredgecolor='k')
    else:
        box(df_, x=col, y=metric, c=create_palette(df_, col, 'Greys'), ax=ax, order=sorted(df_[col].unique()))
    name = col_names_[col] if col in col_names_ else col
    name = name if i in [10,11,12,13,14] else ''
    metric = metric if i in [0,5,10] else ''
    format_ax(ax=ax, xlabel=name, ylabel=metric, reduced_spines=True)

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'key_options_vs_GBC.png'), dpi=500)


##


# More on AF confident detection
metric = 'NMI'
df_agg = df.groupby(['sample', 'af_confident_detection']).apply(lambda x: x[[metric, 'n_cells']].median()).reset_index()# .query('sample=="MDA_clones"')


# Create a figure and a set of subplots with shared x-axis
fig, axs = plt.subplots(1,3,figsize=(12,3.5))

for i,sample in enumerate(['MDA_clones', 'MDA_PT', 'MDA_lung']):

    df_ = df_agg.query('sample==@sample')
    ax = axs.ravel()[i]
    ax.plot(df_['af_confident_detection'], df_[metric], 'bo--', label='NMI')
    ax.set_xlabel('af_confident_detection')
    ax.set_ylabel(metric, color='b')
    ax.tick_params(axis='y', labelcolor='b')
    ax.set(xlabel='Min AF confident detection', title=sample)

    ax_ = ax.twinx()
    ax_.plot(df_['af_confident_detection'], df_['n_cells'], 'gx--', label='n_cells')
    ax_.set_ylabel('n_cells', color='g')
    ax_.tick_params(axis='y', labelcolor='g')

fig.tight_layout()
fig.savefig(os.path.join(path_results, f'AF_{metric}.png'), dpi=500)



##
