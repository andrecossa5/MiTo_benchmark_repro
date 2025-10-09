"""
Fig 5:
- Compare cell phylogenies from different scLT systems
- Visualize:
    - Example cell phylogenies
    - Quantitative metrics to evaluate cell phylogeny inference output
"""

import os
import pandas as pd
import pickle
import mito as mt
import matplotlib
import matplotlib.pyplot as plt
import plotting_utils as plu
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'results', 'others', 'Fig5')
path_figures = os.path.join(path_main, 'results', 'figures', 'Fig5')

# Set visualization params
plu.set_rcParams({'figure.dpi':350})


##


# 1. Choose best combination scLT system / distance metric (within the ones tested) ---------#

# Read data table, format names
df = pd.read_csv(os.path.join(path_data, 'scLT_metrics.csv'), index_col=0)
df = (
    df.drop(columns=['job_id'])
    .pivot_table(index=['method', 'sample'], values='value', columns='metric')
    .reset_index()
)
d = {
    'MAESTER':'weighted_jaccard', 'scWGS':'jaccard',
    'Cas9_unweighted':'hamming', 'Cas9_weighted':'weighted_hamming',
    'RedeeM_unweighted':'jaccard', 'RedeeM_weighted':'weighted_jaccard',
    'EPI-clone':'jaccard'
}
df['metric'] = df['method'].map(d)
df['scLT_system'] = df['method'] 
df.loc[df['method'].str.contains('Cas9'), 'scLT_system'] = 'Cas9'
df.loc[df['method'].str.contains('RedeeM'), 'scLT_system'] = 'RedeeM'

# Select top distance metrics for Cas9 and RedeeM
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


# 2. Fig 5b. Plot 4 representative phylogenies (one per scLT system) ---------#

# Choose phylogenies with richest character matrices
(
    df.groupby('scLT_system')
    .apply(lambda x: x.sort_values('n_characters', ascending=False).head(1))
    [['scLT_system', 'sample']]
)

# Plot
fig, axs = plt.subplots(1,5,figsize=(10.5,2.7))

samples = ['PD34493', '3513_NT_T4', 'LARRY_mouse3', 'Young1.T1.HSP', 'MDA_PT']
method = ['scWGS', 'Cas9', 'EPI-clone', 'RedeeM', 'MAESTER']
d = dict(zip(method, samples))

for i,(method,sample) in enumerate(d.items()):

    path_tree = os.path.join(path_data, 'trees', f'{sample}_tree.pickle')
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
fig.savefig(os.path.join(path_figures, 'Fig_5b.pdf'))


##


# 2. Fig 5c. Visualize scLT comparison results with key evaluation metrics ---------#

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


# Plot
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
fig.savefig(os.path.join(path_figures, 'Fig_5c.pdf'))



##