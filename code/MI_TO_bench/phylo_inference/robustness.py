"""
Do we have robust phylogenies??
"""

import os
from mito_utils.utils import *
from mito_utils.clustering import *
from mito_utils.heatmaps_plots import *
from mito_utils.plotting_base import *
matplotlib.use('macOSX')


##


# Get metrics
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro'
path_data = os.path.join(path_main, 'results', 'MI_TO_bench', 'phylo_inference', 'trees')
path_results = os.path.join(path_main, 'results', 'MI_TO_bench', 'phylo_inference')

# Extract data
# tree_metrics = []
# bench_results = []
# for folder, sub, files in os.walk(path_data):
#     for file in files:
#         splitted = os.path.join(folder, file).split('/')
#         method = splitted[-4]
#         sample = splitted[-3]
#         job_id = splitted[-2]
#         if file.startswith('tree'):
#             tree_metrics.append(pd.read_csv(os.path.join(folder, file), index_col=0).assign(method=method, sample=sample))
#         elif file.startswith('bench'):
#             with open(os.path.join(folder, file), 'rb') as f:
#                 d = pickle.load(f)
#             res = {}
#             res['ARI'] = [ d[k]['ARI'] for k in d ]
#             res['NMI'] = [ d[k]['NMI'] for k in d ]
#             res['NA'] = [ d[k]['NA'] for k in d ]
#             res['sample'] = [ sample for _ in range(len(d)) ]
#             res['job_id'] = [ job_id for _ in range(len(d)) ]
#             res['method'] = list(d.keys())
#             bench_results.append(pd.DataFrame(res))
# df_metrics = (
#     pd.concat(tree_metrics)
#     .pivot_table(values='value', columns='metric', index=['method', 'sample', 'job_id'])
#     .reset_index()
# )
# df_metrics.groupby(['method']).median().T

# Read tables
df_metrics = pd.read_csv(os.path.join(path_data, 'metrics.csv'), index_col=0)
metrics = [
    'median_time', 'n_clones', 'median_n_assigned_char_per_clone', 
    'corr_distances', 'median_CI',
    'median_support', 'median_support_biggest_clades', 'median_support_mut_clades', 
]
df = df_metrics[['method']+metrics]

names = dict(zip(
    metrics, 
    ['Depth', 'n clones', 'n characters per clone', 'Tree vs char- dists correlation',
     'CI', 'Support', 'Support biggest clades', 'Support mutated clades'
    ]
))



fig, axs = plt.subplots(2,4,figsize=(10,5),sharex=True)

order = df.groupby('method')['median_support_biggest_clades'].median().sort_values().index
colors = create_palette(df, 'method', ten_godisnot)

for i,metric in enumerate(metrics):
    ax = axs.ravel()[i]
    box(df, x='method', y=metric, ax=ax, c=colors, order=order)
    format_ax(ax=ax, reduced_spines=True, ylabel=names[metric], rotx=90)

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'mt_phylogeny_robustness.png'), dpi=700)


##
