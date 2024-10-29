"""
Is MI-TO better at solving clones?
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
path_data = os.path.join(path_main, 'results', 'MI_TO_bench', 'benchmark', 'NJ')
path_results = os.path.join(path_main, 'results', 'MI_TO_bench', 'benchmark')

# Extract from bench folder
L = []
for folder, sub, files in os.walk(path_data):
    for file in files:
        if file.startswith('bench'):
            splitted = os.path.join(folder, file).split('/')
            sample = splitted[-3]
            job_id = splitted[-2]
            res = {}
            with open(os.path.join(folder, file), 'rb') as f:
                d = pickle.load(f)
                res['ARI'] = [ d[k]['ARI'] for k in d ]
                res['NMI'] = [ d[k]['NMI'] for k in d ]
                res['NA'] = [ d[k]['NA'] for k in d ]
                res['sample'] = [ sample for _ in range(len(d)) ]
                res['job_id'] = [ job_id for _ in range(len(d)) ]
            L.append(pd.DataFrame(res, index=d.keys()))

df = pd.concat(L).reset_index(names=['method'])


##


# Viz performance
fig, axs = plt.subplots(1,2,figsize=(8,4))

df['sample'] = pd.Categorical(df['sample'], categories=['MDA_clones', 'MDA_PT', 'MDA_lung'])
order = df.groupby('method')['ARI'].median().sort_values().index

colors = { k:v for k,v in zip(order, sc.pl.palettes.vega_10_scanpy) }

box(df, x='sample', y='ARI', by='method', ax=axs[0], c=colors, hue_order=order)
format_ax(ax=axs[0], reduced_spines=True, ylabel='ARI')
box(df, x='sample', y='NMI', by='method', ax=axs[1], c=colors, hue_order=order)
format_ax(ax=axs[1], reduced_spines=True, ylabel='NMI')
add_legend(label='Method', ax=axs[1], colors=colors, loc='center right', 
           bbox_to_anchor=(.55,1.15), ncols=df['method'].unique().size, artists_size=8, label_size=9, ticks_size=8)

fig.subplots_adjust(top=.8, bottom=.1, left=.1, right=.9)
fig.savefig(os.path.join(path_results, 'clonar_reconstruction_performance.png'), dpi=500)


##


# Viz
df['sample'] = df['sample'].astype('str')
df = df.groupby(['sample', 'job_id'])['ARI'].median().sort_values().reset_index()
top_3_jobs = df.groupby('sample').apply(lambda x: x['job_id'].values[x['ARI'].argmax()]).to_dict()


fig, axs = plt.subplots(3,5,figsize=(15,9))

for i,sample in enumerate(['MDA_clones', 'MDA_PT', 'MDA_lung']):

    # Read
    afm = sc.read(os.path.join(path_data, sample, top_3_jobs[sample], 'afm.h5ad'))
    reduce_dimensions(afm, metric='jaccard', bin_method='vanilla')

    with open(os.path.join(path_data, sample, top_3_jobs[sample], 'bench_clonal_recontruction.pickle'), 'rb') as f:
        d = pickle.load(f)
    for method in d:
        labels = d[method]['labels']
        afm.obs[method] = labels
        afm.obs[method][afm.obs[method].isna()] = 'unassigned'

    # Colors
    palette = sc.pl.palettes.vega_20_scanpy if sample != 'MDA_PT' else sc.pl.palettes.godsnot_102
    gbc_colors,_ = assign_matching_colors(afm.obs, 'GBC', 'MI_TO', palette)

    # Plot
    for ax,method in zip(axs[i,:],['GBC', 'MI_TO', 'leiden', 'vireoSNP', 'CClone']):

        if method == 'GBC':
            colors = gbc_colors
        else:
            _,colors = assign_matching_colors(afm.obs, 'GBC', method, palette)
        
        draw_embeddings(
            afm.obs, cat=method, s=10 if sample != 'MDA_clones' else 30, 
            ax=ax, 
            title=f'{method} (n={afm.obs[method].unique().size})',
            legend_kwargs={'colors':colors}, axes_kwargs={'legend':False}
        )


fig.tight_layout()
fig.savefig(os.path.join(path_results, 'viz_clonal_reconstruction.png'), dpi=500)


##
