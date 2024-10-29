"""
Contribute of individual options to individual metrics.
"""

import os
from lightgbm import LGBMRegressor
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
metrics = pd.read_csv(os.path.join(path_data, 'metrics.csv'), header=None).iloc[:,0].to_list()
metrics = [ x for x in metrics if x != 'median_target/untarget_coverage_logratio']
options = pd.read_csv(os.path.join(path_data, 'options.csv'), header=None).iloc[:,0].to_list()
df = df.drop(columns=['median_target/untarget_coverage_logratio', 'cell_subset', 'cell_filter'])


##


# Prep data
df_options = df.loc[:,df.columns.isin(options)].dropna()
variable_options = df_options.nunique().loc[lambda x: x>1].index

X = df_options.loc[:,variable_options]
for x in X.columns:
    try:
        X[x] = X[x].astype(float)
    except:
        X[x] = X[x].astype('category')

##

# Select test metrics
metric_tests = [
    'AUPRC', 'ARI', 'NMI', 'freq_lineage_biased_muts', 
    'n_cells', 'median_n_vars_per_cell', 'density', 'transitions_vs_transversions_ratio'
]
names = {
    'freq_lineage_biased_muts' : "% lineage-biased \n MT-SNVs",
    'n_cells' : 'n cells',
    'median_n_vars_per_cell' : 'n MT-SNVs per cell',
    'density' : 'AFM density',
    'transitions_vs_transversions_ratio' : 'n transitions / \n n transversions'
}


# Here we go
fig, axs = plt.subplots(2,4,figsize=(14,6))

for i, metric in enumerate(metric_tests):

    ax = axs.ravel()[i]
    y = df.loc[df_options.index, metric]

    model = LGBMRegressor(n_estimators=100, learning_rate=0.01, min_data_in_leaf=1, importance_type='split')
    model.fit(X, y)
    
    feature_importances = model.feature_importances_
    df_ = (
        pd.Series(feature_importances, X.columns)
        .sort_values(ascending=False)
        .to_frame('Importance')
    )

    ax.hlines(y=df_.index, xmin=0, xmax=df_['Importance'], color='k', linewidth=1.5, alpha=.7)
    ax.plot(df_['Importance'], df_.index, "o", color='darkred', markersize=8, markeredgecolor='k')
    ax.invert_yaxis()
    title = metric if not metric in names else names[metric]
    format_ax(ax=ax, title=title, reduced_spines=True, xlabel='Feature importance')

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'feature_importance_meta_analysis.png'), dpi=500)


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


# Colors
pp_method_colors = create_palette(df_, 'pp_method', sc.pl.palettes.vega_10)
bin_method_colors = { "vanilla" : '#C9D2CD' , "MI_TO" : '#CA5E00'}


##


fig, ax = plt.subplots(figsize=(4.,4.))
scatter(df_, x='Association with GBC score', y='Yield score', by='pp_method', c=pp_method_colors, s=100,  ax=ax)
format_ax(ax=ax, xlabel='Association with GBC score', ylabel='Yield score')
add_legend(ax=ax, colors=pp_method_colors, loc='upper right', bbox_to_anchor=(1,1), artists_size=8, 
           label_size=9, ticks_size=8, label='Preprocessing\nmethod')
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'overall_scatter.png'), dpi=500)


##


fig, axs = plt.subplots(3,3,figsize=(8,6.5))
df_.columns
metrics_ = [
    'n_cells', 'n_vars', 'median_n_vars_per_cell',
    'n_dbSNP', 'n_REDIdb', 'transitions_vs_transversions_ratio',
    'ARI', 'NMI', 'AUPRC'
]
metric_names_ = {'transitions_vs_transversions_ratio':'Transitions /\nTransversion', 
                 'median_n_vars_per_cell':'Median n_vars per cell'}
for i,metric in enumerate(metrics_):
    ax = axs.ravel()[i]
    order = df_.groupby('pp_method')[metric].median().sort_values().index
    x = df_.groupby('pp_method')[metric].median().sort_values()
    df_feat = df.groupby(groupings, dropna=False)[metric].median().reset_index()
    name = metric_names_[metric] if metric in metric_names_ else metric
    box(df_feat, x='pp_method', y=metric, ax=ax, c=pp_method_colors, order=order)
    format_ax(ax=ax,  xticks=[], title=name, reduced_spines=True)

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'metrics_by_method.png'), dpi=500)



##


fig, ax = plt.subplots(figsize=(4.,4.))
scatter(df_, x='Association with GBC score', y='Yield score', by='bin_method', c=bin_method_colors, s=100,  ax=ax)
format_ax(ax=ax, xlabel='Association with GBC score', ylabel='Yield score')
add_legend(ax=ax, colors=bin_method_colors, loc='upper right', bbox_to_anchor=(1,1), artists_size=8, 
           label_size=9, ticks_size=8, label='Binarization\nmethod')
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'bin_scatter.png'), dpi=500)


##


fig, axs = plt.subplots(3,3,figsize=(8,6.5), sharex=True)
df_.columns
metrics_ = [
    'n_cells', 'n_vars', 'median_n_vars_per_cell',
    'n_dbSNP', 'n_REDIdb', 'transitions_vs_transversions_ratio',
    'ARI', 'NMI', 'AUPRC'
]
metric_names_ = {'transitions_vs_transversions_ratio':'Transitions /\nTransversion', 
                 'median_n_vars_per_cell':'Median n_vars per cell'}

order = df_.groupby('pp_method')['ARI'].mean().sort_values().index

for i,metric in enumerate(metrics_):
    ax = axs.ravel()[i]
    df_feat = df.groupby(groupings, dropna=False)[metric].mean().reset_index()
    name = metric_names_[metric] if metric in metric_names_ else metric
    box(df_feat, x='pp_method', y=metric, ax=ax, by='bin_method', c=bin_method_colors, hue_order=bin_method_colors.keys(), saturation=1, order=order)
    format_ax(ax=ax,  title=name, reduced_spines=True, rotx=90)

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'metrics_by_bin.png'), dpi=500)


##


