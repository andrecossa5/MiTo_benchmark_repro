"""
Supp Fig ... distance metrics comparison.
Impact of different distance metrics on several metrics.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import mito as mt
import matplotlib
import matplotlib.pyplot as plt
import plotting_utils as plu
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'data', 'tune')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')
path_results = os.path.join(path_main, 'results', 'others', 'Supp')


##


# 2. Single jobs. Only maegatk, filter2, MiTo. -------------------------- 
path_data = os.path.join(path_main, 'results', 'others', 'Fig2')

# Read selected jobs
L = []
for x in os.listdir(path_data):
    if x.endswith('filtered_jobs.csv'):
        L.append(pd.read_csv(os.path.join(path_data, x), index_col=0))
df = pd.concat(L)

# Filter only: pp_method == maegatk, cell_filter==filter2, filtering==MiTo
df = (
    df.loc[
        (df['cell_filter']=='filter2') & \
        (df['pp_method'].isin(['maegatk'])) & \
        (df['filtering']=='MiTo') & \
        (df['min_cell_number']==10)
    ]
)

##

df.groupby('metric')[['ARI', 'NMI', 'corr']].describe().T



fig, axs = plt.subplots(1,3,figsize=(6,3.5),sharey=True)

metric = 'AUPRC'
bin_method_colors = {
    'jaccard': 'orange',
    'weighted_jaccard': 'green'
}

for i,sample in enumerate(['MDA_clones', 'MDA_PT', 'MDA_lung']):
    ax = axs[i]
    plu.box(df.query('sample==@sample'), 'af_confident_detection', metric, 
            categorical_cmap=bin_method_colors, by='metric', by_order=['jaccard', 'weighted_jaccard'], ax=ax)#, alpha=1)
    plu.format_ax(ax=ax, title=sample, reduced_spines=True, xlabel='Min confident AF', ylabel=metric)

fig.tight_layout()
plt.show()
fig.savefig(os.path.join(path_figures, f'Supp_fig_binarization_{metric}_selected_jobs.pdf'))


##
