"""
Visualize MT- and GBC- clone dynamics.
"""

import os
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import mito as mt
import plotting_utils as plu
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import norm
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'results', 'others', 'Fig4', 'longitudinal')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')

# Read data
afm = sc.read(os.path.join(path_data, 'afm_filtered.h5ad'))
tree_metrics = pd.read_csv(os.path.join(path_data, 'tree_metrics.csv'), index_col=0)
with open(os.path.join(path_data, 'annotated_tree.pickle'), 'rb') as f:
    tree = pickle.load(f)

##


def gaussian_smooth(x, y, grid, sd):
    weights = np.transpose([norm.pdf(grid, m, sd) for m in x])
    weights = weights / weights.sum(0)
    return (weights * y).sum(1)


##


# GBC clones
clone_type = 'GBC'
colors = plu.create_palette(tree.cell_meta, clone_type, sc.pl.palettes.default_102)

df_long = (
    tree.cell_meta
    [['MiTo clone', 'GBC', 'sample']]
    .loc[lambda x: ~x['MiTo clone'].isna()]
    .groupby('sample')
    .apply(lambda x: x[clone_type].value_counts(normalize=True))
    .reset_index()
    .pivot_table(columns='sample', index=clone_type, values='proportion', fill_value=0)
    #.sort_values('MDA_PT', ascending=False)
)

# Viz
plu.set_rcParams()

x = np.arange(df_long.shape[1])
y = [ np.array(x) for x in df_long.values.tolist() ]
colors_ = [ colors[x] for x in df_long.index ] 
grid = np.linspace(-1.3, 4, num=500)
y_smoothed = [ gaussian_smooth(x, y_, grid, .35) for y_ in y ]

##

fig, ax = plt.subplots(figsize=(4.5,2))
ax.stackplot(grid, y_smoothed, baseline="sym", colors=colors_)
plu.format_ax(ax, xticks=['PT', 'lung'], yticks=[])
ax.spines[['right', 'bottom', 'top', 'left']].set_visible(False)
for l in x:    
    ax.axvline(x=l, color='k', linewidth=.5, linestyle='dashed')
fig.tight_layout()
fig.savefig(os.path.join(path_figures, f'{clone_type}_fishplot_long.pdf'))

##