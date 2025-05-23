"""
Plot scLT_systems comparisons.
"""

import os
import pickle
import pandas as pd
import mito as mt
import plotting_utils as plu
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('macOSX')


# Paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro/'
path_results = os.path.join(path_main, 'results/others/Fig5')
path_figures = os.path.join(path_main, 'results/figures/Fig5')

# Read data table, format names
df = pd.read_csv(os.path.join(path_results, 'scLT_metrics.csv'), index_col=0)
df = (
    df.drop(columns=['job_id'])
    .pivot_table(index=['method', 'sample'], values='value', columns='metric')
    .reset_index()
)
d = {
    'MAESTER':'weighted_jaccard', 'scWGS':'jaccard',
    'Cas9_unweighted':'hamming', 'Cas9_weighted':'weighted_hamming',
    'RedeeM_unweighted':'jaccard', 'RedeeM_weighted':'weighted_jaccard',
}
df['metric'] = df['method'].map(d)
df['scLT_system'] = df['method'] 
df.loc[df['method'].str.contains('Cas9'), 'scLT_system'] = 'Cas9'
df.loc[df['method'].str.contains('RedeeM'), 'scLT_system'] = 'RedeeM'

# Select top for Cas9 and RedeeM
(
    df.query('scLT_system=="Cas9" or scLT_system=="RedeeM"')
    .groupby(['scLT_system', 'metric'])
    [['corr_distances', 'median_support', 'median_support_biggest_clades']]
    .median()
)
df = ( 
    df.query('method!="RedeeM_weighted" and method!="Cas9_unweighted"')
    .drop(columns=['method'])
)


##


# Select metrics to plot
metrics = [
    'n_cells', 'n_characters', 'median_char_per_cell', 'density',
    'corr_distances', 'median_support', 'median_support_biggest_clades',
    'frac_remaining_nodes_after_mutationless_edges_collapse'
]
df = df[['scLT_system', 'sample']+metrics]
names = dict(zip(
    metrics, 
    ['n cells', 'n characters', 'n characters per cell', 'Density',
     'Tree vs char- dists correlation',
     'Support', 'Support biggest clades', '% nodes after tree pruning'
    ]
))


##


# Viz
plu.set_rcParams({'figure.dpi':300})

fig, axs = plt.subplots(2,4,figsize=(8.5,4.5),sharex=True)

order = df.groupby('scLT_system')['median_support_biggest_clades'].median().sort_values().index
colors = plu.create_palette(df, 'scLT_system', 'Spectral_r')

for i,metric in enumerate(metrics):
    ax = axs.ravel()[i]
    plu.violin(df, x='scLT_system', y=metric, ax=ax, categorical_cmap=colors, x_order=order, 
               kwargs={'alpha':.8, 'bw_adjust':.6})
    plu.strip(df, x='scLT_system', y=metric, ax=ax, categorical_cmap=colors, x_order=order)
    plu.format_ax(ax=ax, reduced_spines=True, ylabel=names[metric], xlabel='', rotx=90)

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'scLT_system_overall_comparison.pdf'))


##


# Choose top samples
(
    df.groupby('scLT_system')
    .apply(lambda x: x.sort_values('n_characters', ascending=False).head(1))
    [['scLT_system', 'sample']]
)


# Viz trees
plu.set_rcParams({'figure.dpi':150})

fig, axs = plt.subplots(1,4,figsize=(8.5,2.7))

samples = ['PD34493', '3513_NT_T4', 'Young1.T1.HSP', 'MDA_PT']
method = ['scWGS', 'Cas9', 'RedeeM', 'MAESTER']
d = dict(zip(method, samples))

for i,(method,sample) in enumerate(d.items()):
    path_tree = os.path.join(path_results, 'trees', f'{sample}_tree.pickle')
    with open(path_tree, 'rb') as f:
        tree = pickle.load(f)

    n_cells, n_vars = tree.layers['transformed'].shape
    if method == 'Cas9':
        n_vars = tree.layers['transformed'].nunique().sum() - \
                 tree.layers['transformed'].shape[1]

    ax = axs.ravel()[i]
    mt.pl.plot_tree(tree, ax=ax, orient=90, branch_kwargs={'linewidth':.5})
    ax.set(title=f'{method}: {sample}\nn leaves: {n_cells}\nn characters: {n_vars}')

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'trees.pdf'))


##