"""
Binomial mixture model visualization
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import plotting_utils as plu
from mito.ut.stats_utils import get_components
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
sample = 'MDA_clones'
job = '77feb9f505'
path_data = os.path.join(path_main, 'data', 'lineage_inference', 'UPMGA', sample, job)
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')


##


# Load data
afm = sc.read(os.path.join(path_data, 'afm_filtered.h5ad'))
AD = afm.layers['AD'].toarray()
DP = afm.layers['site_coverage'].toarray()

# Choose idx
# afm.var['n5']
idx = 4

# Viz single MT fitted components
ad = AD[:,idx]
dp = DP[:,idx]
x0, x1 = get_components(ad, dp)
# Re=factor as df
df = pd.DataFrame(
    { 
        'value' : np.concatenate([ad, x0, x1]), 
        'type' : ['Empirical']*ad.size + ['C0']*x0.size + ['C1']*x1.size
    }
)

# Plot
plu.set_rcParams({'figure.dpi':80})
fig, ax = plt.subplots(figsize=(3.5,3))
cmap = {'Empirical' : 'grey', 'C0' : '#D58E8E', 'C1' : '#A13A3A'}
plu.dist(df, x='value', by='type', ax=ax, categorical_cmap=cmap, alpha=.2, fill=True, linewidth=1, bw_adjust=3)
plu.format_ax(ax=ax, xlabel='nUMIs ALT allele', ylabel='Density', reduced_spines=True)
plu.add_legend(ax=ax, colors=cmap, label='', loc='upper right', bbox_to_anchor=(1,1))
plt.tight_layout()
fig.savefig(os.path.join(path_figures, 'binomial_mixture_example.pdf'), dpi=300)


##
