"""
Is MiTo better at solving clones than other clonal reconstruction approaches?
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
path_data = os.path.join(path_main, 'results', 'MI_TO_bench', 'phylo_inference_final', 'refining_mito_maegatk')
path_results = os.path.join(path_main, 'results', 'MI_TO_bench', 'phylo_inference_final')
path_colors = os.path.join(path_main, 'data/MI_TO_bench/clones_colors_sc.pickle')

# Extract data
tree_metrics = []
bench_results = []
for folder, sub, files in os.walk(path_data):
    for file in files:
        splitted = os.path.join(folder, file).split('/')
        # method = splitted[-4]
        sample = splitted[-3]
        job_id = splitted[-2]
        if file.startswith('tree'):
            tree_metrics.append(pd.read_csv(os.path.join(folder, file), index_col=0).assign(sample=sample)) # method=method
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


##


# MiTo metrics           
df_metrics = (
    pd.concat(tree_metrics)
    .pivot_table(values='value', columns='metric', index=['sample', 'job_id'])
    .reset_index()
)
df_metrics.groupby('sample').median().T

# MiTo vs others bench
df_bench = pd.concat(bench_results)
df_bench.groupby(['sample', 'method']).median()


##


# Viz performance
fig, axs = plt.subplots(1,2,figsize=(9,4))

df_bench['sample'] = pd.Categorical(df_bench['sample'], categories=['MDA_clones', 'MDA_lung', 'MDA_PT'])
order = df_bench.groupby('method')['ARI'].median().sort_values().index

colors = { k:v for k,v in zip(order, sc.pl.palettes.vega_10_scanpy) }

sns.stripplot(data=df_bench, x='sample', y='ARI', palette=colors.values(), 
              hue='method', hue_order=order, size=5, ax=axs[0], dodge=True, edgecolor='k', linewidth=.5)
sns.barplot(data=df_bench, x='sample', y='ARI', palette=colors.values(), 
            hue='method', hue_order=order, ax=axs[0], dodge=True)
axs[0].get_legend().remove()
axs[0].set_ylim((-.01,1))
axs[0].axhline(.9, linestyle='--', color='k')
format_ax(ax=axs[0], xlabel='', ylabel='ARI', reduced_spines=True, xlabel_size=10, ylabel_size=10)

sns.stripplot(data=df_bench, x='sample', y='NMI', palette=colors.values(), 
              hue='method', hue_order=order, size=5, ax=axs[1], dodge=True, edgecolor='k', linewidth=.5)
format_ax(ax=axs[1], xlabel='', ylabel='NMI', reduced_spines=True, xlabel_size=10, ylabel_size=10)
sns.barplot(data=df_bench, x='sample', y='NMI', palette=colors.values(), 
            hue='method', hue_order=order, ax=axs[1], dodge=True)
axs[1].get_legend().remove()
axs[1].set_ylim((-.01,1))
axs[1].axhline(.9, linestyle='--', color='k')
format_ax(ax=axs[1], xlabel='', ylabel='NMI', reduced_spines=True, xlabel_size=10, ylabel_size=10)
add_legend(ax=axs[1], colors=colors, loc='upper left', bbox_to_anchor=(1,1), artists_size=10, label='Method', label_size=10, ticks_size=10)

fig.subplots_adjust(right=.75, top=.85, left=.15, bottom=.15)
fig.savefig(os.path.join(path_results, 'clonal_reconstruction_performance.png'), dpi=700)


##


# Viz
df_bench['sample'] = df_bench['sample'].astype('str')
df_bench = df_bench.groupby(['sample', 'job_id'])['ARI'].median().sort_values().reset_index()
top_3_jobs = df_bench.groupby('sample').apply(lambda x: x['job_id'].values[x['ARI'].argmax()]).to_dict()

# Colors
with open(path_colors, 'rb') as f:
    clone_colors = pickle.load(f)

# Viz
fig, axs = plt.subplots(3,5,figsize=(14,9))

for i,sample in enumerate(['MDA_clones', 'MDA_PT', 'MDA_lung']):

    # Read
    afm = sc.read(os.path.join(path_data, sample, top_3_jobs[sample], 'afm.h5ad'))
    reduce_dimensions(afm)

    with open(os.path.join(path_data, sample, top_3_jobs[sample], 'bench_clonal_recontruction.pickle'), 'rb') as f:
        d = pickle.load(f)
    
    for method in ['MiTo', 'vireoSNP', 'leiden', 'CClone']:
        labels = d[method]['labels']
        afm.obs[method] = labels
        afm.obs[method][afm.obs[method].isna()] = 'unassigned'
        afm.obs[method] = pd.Categorical(labels)

    # Colors
    n_gt = afm.obs['GBC'].cat.categories.size
    palette = list(clone_colors.values())
    gbc_colors,_ = assign_matching_colors(afm.obs, 'GBC', 'MiTo', palette)
    afm.uns['GBC_colors'] = list(gbc_colors.values())

    # Plot
    for ax,method in zip(axs[i,:],['GBC', 'MiTo', 'vireoSNP', 'leiden', 'CClone']):

        if method == 'GBC':
            colors = gbc_colors
            title = f'{method}\nn clones={n_gt}'
        else:
            ari = d[method]['ARI']
            n_clones = afm.obs[method].unique().size
            _,colors = assign_matching_colors(afm.obs, 'GBC', method, palette)
            title = f'{method}\nn clones={n_clones}, ARI={ari:.2f}'
        
        afm.uns[f'{method}_colors'] = list(colors.values())
        sc.pl.embedding(
            afm, basis='X_umap', color=method, ax=ax,
            frameon=False, save=False, legend_loc=None, show=False, 
            size=40 if sample != 'MDA_clones' else 75,
            title=title
        )

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'viz_clonal_reconstruction.png'), dpi=500)


##