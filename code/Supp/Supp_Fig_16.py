"""
Supp Fig 16:
- Read PT-lung dataset filtered AFM
- Visualize MT-SNVs and cell-cell distances with clustered heatmaps.
"""

import os
import pickle
import scanpy as sc
import mito as mt
import plotting_utils as plu
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_lineage_inference = os.path.join(path_main, 'results', 'others', 'Fig4', 'lineage_inference')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')


# Set visualization params
plu.set_rcParams({'figure.dpi':350})


##


# Fig 16 a,b. AFM and cell-cell distances visualization ------------------------------------#

# Read filtered AFM cell phylogeny
afm = sc.read(os.path.join(path_lineage_inference, 'afm_filtered.h5ad'))
with open(os.path.join(path_lineage_inference, 'annotated_tree.pickle'), 'rb') as f:
    tree = pickle.load(f)


##


# Distances
fig, ax = plt.subplots(figsize=(4,4))
mt.pl.heatmap_distances(afm, tree=tree, ax=ax)
fig.tight_layout()
plu.save_best_pdf_quality(fig, figsize=(4,4), path=path_figures, name='Supp_Fig_16a.pdf', dpi=1000)

# Variants
fig, ax = plt.subplots(figsize=(4,4))
mt.pl.heatmap_variants(afm, tree=tree, ax=ax, kwargs={'x_names_size':3}, cmap='afmhot_r')
fig.tight_layout()
plu.save_best_pdf_quality(fig, figsize=(4,4), path=path_figures, name='Supp_Fig_16b.pdf', dpi=1000)


##