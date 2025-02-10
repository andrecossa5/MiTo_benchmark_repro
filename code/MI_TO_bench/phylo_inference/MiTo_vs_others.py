"""
Is MI-TO better at solving clones than other clonal reconstruction approaches?
"""

import os
from mito_utils.utils import *
from mito_utils.plotting_base import *
from mito_utils.dimred import *
from mito_utils.embeddings_plots import *
matplotlib.use('macOSX')


##


# Get metrics
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro'
path_data = os.path.join(path_main, 'results', 'MI_TO_bench', 'phylo_inference', 'NJ')
path_results = os.path.join(path_main, 'results', 'MI_TO_bench', 'phylo_inference')
path_colors = os.path.join(path_main, 'data/MI_TO_bench/clones_colors_sc.pickle')



len(list(os.walk(path_data)))


# Extract data
tree_metrics = []
bench_results = []
for folder, sub, files in os.walk(path_data):
    for file in files:
        splitted = os.path.join(folder, file).split('/')
        method = splitted[-4]
        sample = splitted[-3]
        job_id = splitted[-2]
        if file.startswith('tree'):
            tree_metrics.append(pd.read_csv(os.path.join(folder, file), index_col=0).assign(method=method, sample=sample))
        elif file.startswith('bench'):
            with open(os.path.join(folder, file), 'rb') as f:
                d = pickle.load(f)
            res = {}
            res['ARI'] = [ d[k]['ARI'] for k in d ]
            res['NMI'] = [ d[k]['NMI'] for k in d ]
            res['% unassigned'] = [ d[k]['% unassigned'] for k in d ]
            res['sample'] = [ sample for _ in range(len(d)) ]
            res['job_id'] = [ job_id for _ in range(len(d)) ]
            res['method'] = list(d.keys())
            bench_results.append(pd.DataFrame(res))
df_metrics = (
    pd.concat(tree_metrics)
    .pivot_table(values='value', columns='metric', index=['method', 'sample', 'job_id'])
    .reset_index()
)
df_metrics.groupby(['method']).median().T
df_bench = pd.concat(bench_results)

df_bench.groupby(['sample', 'method']).mean()

# Read tables
df_metrics = pd.read_csv(os.path.join(path_data, 'metrics.csv'), index_col=0)
df_metrics['method'] = df_metrics['method'].map(lambda x: f'MiTo-{x}')
df_bench = pd.read_csv(os.path.join(path_data, 'bench_df.csv'), index_col=0)
cols = ['sample', 'job_id', 'ARI', 'NMI', 'method']
df_bench = pd.concat([ df_bench.query('method!="MI_TO"')[cols], df_metrics[cols] ])


##


# Viz performance
fig, axs = plt.subplots(1,2,figsize=(10,4))

df_bench['sample'] = pd.Categorical(df_bench['sample'], categories=['MDA_clones', 'MDA_lung', 'MDA_PT'])
order = df_bench.groupby('method')['ARI'].median().sort_values().index

colors = { k:v for k,v in zip(order, ten_godisnot) }

box(df_bench, x='sample', y='ARI', by='method', ax=axs[0], c=colors, hue_order=order)
format_ax(ax=axs[0], reduced_spines=True, ylabel='ARI')
box(df_bench, x='sample', y='NMI', by='method', ax=axs[1], c=colors, hue_order=order)
format_ax(ax=axs[1], reduced_spines=True, ylabel='NMI')
add_legend(label='Method', ax=axs[1], colors=colors, loc='upper left', 
           bbox_to_anchor=(1,1), ncols=1,#round(df_bench['method'].unique().size / 2), 
           artists_size=8, label_size=9, ticks_size=8)

fig.subplots_adjust(top=.9, bottom=.1, left=.1, right=.75)
# fig.tight_layout()
plt.show()
fig.savefig(os.path.join(path_results, '../clonal_reconstruction_performance.png'), dpi=700)


##


# Viz
df_bench['sample'] = df_bench['sample'].astype('str')
df_bench = df_bench.groupby(['sample', 'job_id'])['ARI'].median().sort_values().reset_index()
top_3_jobs = df_bench.groupby('sample').apply(lambda x: x['job_id'].values[x['ARI'].argmax()]).to_dict()
# top_3_jobs

# Colors
with open(path_colors, 'rb') as f:
    clone_colors = pickle.load(f)

# Viz
fig, axs = plt.subplots(1,5,figsize=(12,3))

for i,sample in enumerate(['MDA_clones']):#  'MDA_PT', 'MDA_lung']):

    # Read
    afm = sc.read(os.path.join(path_data, 'top_spaces', top_3_jobs[sample], 'afm.h5ad'))
    reduce_dimensions(afm)

    with open(os.path.join(path_data, 'top_spaces', top_3_jobs[sample], 'bench_clonal_recontruction.pickle'), 'rb') as f:
        d = pickle.load(f)
    for method in d:
        labels = d[method]['labels']
        method = method if method != 'MI_TO' else 'MiTo-UPMGA'
        afm.obs[method] = labels
        afm.obs[method][afm.obs[method].isna()] = 'unassigned'

    # Colors
    palette = list(clone_colors.values())
    gbc_colors,_ = assign_matching_colors(afm.obs, 'GBC', 'MiTo-UPMGA', palette)

    # Plot
    for ax,method in zip(axs,['GBC', 'MiTo-UPMGA', 'vireoSNP', 'leiden', 'CClone']):
        if method == 'GBC':
            colors = gbc_colors
        else:
            _,colors = assign_matching_colors(afm.obs, 'GBC', method, palette)
        draw_embeddings(
            afm.obs, cat=method, s=13 if sample != 'MDA_clones' else 40, 
            ax=ax, 
            title=f'{method} (n={afm.obs[method].unique().size})',
            legend_kwargs={'colors':colors}, axes_kwargs={'legend':False}
        )

fig.tight_layout()
plt.show()
# fig.savefig(os.path.join(path_results, '../viz_clonal_reconstruction.png'), dpi=500)


##


# Improves over current metrics

# Read data
from mito_utils.metrics import *


sample = 'MDA_PT'
afm = sc.read(os.path.join(path_data, 'top_spaces', top_3_jobs[sample], 'afm.h5ad'))
reduce_dimensions(afm)
with open(os.path.join(path_data, 'top_spaces', top_3_jobs[sample], 'bench_clonal_recontruction.pickle'), 'rb') as f:
    d = pickle.load(f)

# Get labels
for method in d:
    res = d[method]
    afm.obs[method] = res['labels']

df_labels = afm.obs[['GBC', 'leiden', 'vireoSNP','MI_TO', 'CClone']]

np.sum(df_labels['MI_TO']=='Undefined')

custom_ARI(df_labels.query('MI_TO != "Undefined"')['GBC'], df_labels.query('MI_TO != "Undefined"')['MI_TO'])
custom_ARI(df_labels.loc[lambda x: ~x['vireoSNP'].isna()]['GBC'], df_labels.loc[lambda x: ~x['vireoSNP'].isna()]['vireoSNP'])
normalized_mutual_info_score(df_labels.query('MI_TO != "Undefined"')['GBC'], df_labels.query('MI_TO != "Undefined"')['MI_TO'])
normalized_mutual_info_score(df_labels.loc[lambda x: ~x['vireoSNP'].isna()]['GBC'], df_labels.loc[lambda x: ~x['vireoSNP'].isna()]['vireoSNP'])


##
