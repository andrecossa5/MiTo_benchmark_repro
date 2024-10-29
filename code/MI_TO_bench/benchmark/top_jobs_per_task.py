"""
Contribute of individual options to individual metrics.
"""

import os
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

# Overall
weights = {
    'Mutation Quality': 0.1,
    'Association with GBC': 0.5,
    'Noise robustness' : .2,
    'Connectivity' : .0,
    'Variation' : .0,
    'Yield' : .2
}
metric_annot = {
    'Mutation Quality' : ['n_dbSNP', 'n_REDIdb', 'transitions_vs_transversions_ratio'],
    'Association with GBC' : ['freq_lineage_biased_muts',  'AUPRC', 'ARI', 'NMI'],                               
    'Noise robustness' : ['corr'],
    'Connectivity' : ['density', 'transitivity', 'average_path_length', 'average_degree', 'proportion_largest_component'],
    'Variation' : ['genomes_redundancy', 'median_n_vars_per_cell'],                                                           
    'Yield' : ['n_GBC_groups', 'n_cells']                                                                
} 
groupings = ['job_id']


##


# Relevant metrics
metrics_ = [
    'n_cells', 'n_vars', 'median_n_vars_per_cell',
    'ARI', 'NMI', 'AUPRC'
]
metric_names_ = {'median_n_vars_per_cell':'Median n_vars per cell'}

L = []
samples = ["MDA_clones", 'MDA_lung', 'MDA_PT']
for sample in samples:
    df_ = rank_items(df.query('sample==@sample'), groupings, metrics, weights, metric_annot)
    L.append(df_[metrics_].drop_duplicates().head(50).assign(sample=sample))

df_ = pd.concat(L)


# Annotate
fig, axs = plt.subplots(2,3,figsize=(8,4.5))

sample_colors = create_palette(df_, 'sample', "Reds")

for i,metric in enumerate(metrics_):
    ax = axs.ravel()[i]
    name = metric_names_[metric] if metric in metric_names_ else metric
    xticks = samples if i>2 else []
    box(df_, x='sample', y=metric, ax=ax, c=sample_colors, saturation=.7, order=samples)
    format_ax(ax=ax, title=name, xticks=xticks, reduced_spines=True, rotx=90)
    if i>2:
        ax.set(ylim=(0,1))

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'top_models_by_sample.png'), dpi=500)


##