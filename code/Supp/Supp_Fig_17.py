"""
Supp Fig 17:
- Fishplot MT- and GBC-clones in longitudinal dataset.
"""

import os
import pickle
import numpy as np
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
path_lineage_inference = os.path.join(path_main, 'results', 'others', 'Fig4', 'lineage_inference')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')


# Set visualization params
plu.set_rcParams({'figure.dpi':350})


##


# Utils
def gaussian_smooth(x, y, grid, sd):
    weights = np.transpose([norm.pdf(grid, m, sd) for m in x])
    weights = weights / weights.sum(0)
    return (weights * y).sum(1)


##


# Fig 17. AFM and cell-cell distances visualization ------------------------------------#

# Read filtered AFM cell phylogeny
afm = sc.read(os.path.join(path_lineage_inference, 'afm_filtered.h5ad'))
with open(os.path.join(path_lineage_inference, 'annotated_tree.pickle'), 'rb') as f:
    tree = pickle.load(f)



##


# Prep data
clone_type = 'MiTo clone'  # Change MiTo clone for MT clones
df_long = (
    tree.cell_meta
    [['MiTo clone', 'GBC', 'sample']]
    .loc[lambda x: ~x['MiTo clone'].isna()]
    .groupby('sample')
    .apply(lambda x: x[clone_type].value_counts(normalize=True))
    .reset_index()
    .pivot_table(columns='sample', index=clone_type, values='proportion', fill_value=0)
)

# Colors
colors = plu.create_palette(tree.cell_meta, clone_type, sc.pl.palettes.default_102)
x = np.arange(df_long.shape[1])
y = [ np.array(x) for x in df_long.values.tolist() ]
colors_ = [ colors[x] for x in df_long.index ] 
grid = np.linspace(-1.3, 4, num=500)
y_smoothed = [ gaussian_smooth(x, y_, grid, .35) for y_ in y ]


##


# Fishplot
fig, ax = plt.subplots(figsize=(4.5,2))
ax.stackplot(grid, y_smoothed, baseline="sym", colors=colors_)
plu.format_ax(ax, xticks=['PT', 'lung'], yticks=[])
ax.spines[['right', 'bottom', 'top', 'left']].set_visible(False)
for l in x:    
    ax.axvline(x=l, color='k', linewidth=.5, linestyle='dashed')
fig.tight_layout()
fig.savefig(os.path.join(path_figures, f'Supp_Fig_17_{clone_type}.pdf'))


##