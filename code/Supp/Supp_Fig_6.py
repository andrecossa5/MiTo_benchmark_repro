"""
Supp Fig 6
MiTo (raw sequencing reads pre-processing) vs maegatk (coverage across MT genome).
"""

import os
import pandas as pd
import mito as mt
import matplotlib
import matplotlib.pyplot as plt
import plotting_utils as plu
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'data', 'bench', 'MiTo_vs_maegatk_consensus')
path_QC_cells = os.path.join(path_main, 'data', 'general', 'QC', 'MDA_clones.txt')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')


# Set visualization params
plu.set_rcParams({'figure.dpi':350})



##


# Gather coverage info
maegatk_mock = pd.read_csv(os.path.join(path_data, 'mock_maegatk', 'coverage.txt.gz'), header=None)
maegatk_mock.columns = ['pos', 'cell', 'n']
mito_prep = pd.read_csv(os.path.join(path_data, 'mito_preprocessing', 'coverage.txt.gz'), header=None)
mito_prep.columns = ['pos', 'cell', 'n']
maegatk_mock_wide = maegatk_mock.pivot_table(index='cell', columns='pos', values='n').fillna(0)
mito_prep_wide = mito_prep.pivot_table(index='cell', columns='pos', values='n').fillna(0)
sites = list(set(mito_prep_wide.columns) & set(maegatk_mock_wide.columns))
maegatk_mock_wide = maegatk_mock_wide[sites]
mito_prep_wide = mito_prep_wide[sites]

# Gather qualified cells and subset both coverage table for them
cells = pd.read_csv(path_QC_cells, header=None)[0].unique().tolist()
cells = [ x.split('_')[0] for x in  cells ] 
subset = list(set(cells) & set(maegatk_mock_wide.index) & set(mito_prep_wide.index))

# Plot
fig, axs = plt.subplots(1,2,figsize=(10,5), subplot_kw={'projection': 'polar'})
mt.pl.MT_coverage_by_gene_polar(mito_prep.query('cell in @subset'), sample='MiTo', ax=axs[0])
mt.pl.MT_coverage_by_gene_polar(maegatk_mock.query('cell in @subset'), sample='maegatk', ax=axs[1])
fig.subplots_adjust(top=.8, bottom=.2, left=.2, right=.8 )
fig.savefig(os.path.join(path_figures, 'Supp_Fig_6.pdf'))


##