"""
Supp Fig 10
Supervised characterization of ground-truth clones MT-SNVs load.
"""

import os
import numpy as np
import pandas as pd
import mito as mt
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import plotting_utils as plu
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'data', 'general', 'AFMs')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')


# Set visualization params
plu.set_rcParams({'figure.dpi':350})


##


# 1. Gather n of MT-SNVs per lentiviral clone (GBC labels) -----------------------------------#

L = []
fdr_treshold = .05

samples = ['MDA_clones', 'MDA_PT', 'MDA_lung']
for sample in samples:
    # Read unfiltered AFM
    afm = sc.read(os.path.join(path_data, sample, 'afm_unfiltered.h5ad'))
    afm = mt.pp.filter_cells(afm, cell_filter='filter2')
    # Filter GBC-enriched muts
    afm = mt.pp.filter_afm(afm, lineage_column='GBC', min_cell_number=5, compute_enrichment=True)
    fdrs = afm.var.loc[:,afm.var.columns.str.startswith('FDR')]
    df = pd.DataFrame({
        'cat' : np.where(np.sum(fdrs<=fdr_treshold)>0, 'Defined', 'Undefined'),
        'sample' : [ sample for _ in range(fdrs.columns.size) ],
        'n_muts' : np.sum(fdrs<=fdr_treshold)
    })
    L.append(df)

df_plot = pd.concat(L).reset_index(drop=True)
df_plot['sample'] = pd.Categorical(df_plot['sample'], categories=samples[::-1])


##


# 1. Plot -----------------------------------#

fig, ax = plt.subplots(figsize=(4,3))

cmap = plu.create_palette(df_plot, 'cat', order=['Undefined', 'Defined'], palette='Greys')
plu.bb_plot(df_plot, cov1='sample', cov2='cat', ax=ax, categorical_cmap=cmap, show_y=True)
plu.add_legend(cmap, label='Clone type', loc='center', bbox_to_anchor=(.5,1.4), ncols=2, ax=ax)
fig.subplots_adjust(top=.6, bottom=.3, right=.85, left=.25)
fig.savefig(os.path.join(path_figures, 'Supp_Fig_10a.pdf'))


##


fig, ax = plt.subplots(figsize=(4,3))

cmap = plu.create_palette(df_plot, 'sample', order=samples, palette='tab10')
plu.dist(df_plot, x='n_muts', by='sample', ax=ax, linewidth=.75, alpha=.3, categorical_cmap=cmap)
medians = {}
for sample in samples:
    median = np.median(df_plot.query('sample==@sample')['n_muts'])
    ax.axvline(x=median, linewidth=.75, linestyle='--', c=cmap[sample])
    medians[sample] = median

ax.text(.43, .75, str(int(medians['MDA_clones'])), transform=ax.transAxes, c=cmap['MDA_clones'])
ax.text(.38, .85, str(int(medians['MDA_lung'])), transform=ax.transAxes, c=cmap['MDA_lung'])
ax.text(.33, .95, str(int(medians['MDA_PT'])), transform=ax.transAxes, c=cmap['MDA_PT'])

plu.format_ax(ax=ax, xlabel='n enriched MT-SNVs', ylabel='Density', 
              title=f'n lentiviral clones: {df_plot.shape[0]}', reduced_spines=True)
plu.add_legend(cmap, label='Sample', loc='upper right', bbox_to_anchor=(1,1), ax=ax)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'Supp_Fig_10b.pdf'))


##