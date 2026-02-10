"""
Supp Fig 20a
- Re-analysis of Weng et. al, 2024 dataset: mouse batch 1.
- Visualize n of Cas9 clones, cells and SH.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import plotting_utils as plu
matplotlib.use('macOSX')


# Set visualization params
plu.set_rcParams({'figure.dpi':120})


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'data', 'redeem_mouse')
path_results = os.path.join(path_main, 'results', 'others', 'Supp20')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')


# Gather cells meta
samples = ['DKO_1', 'CD47_1', 'Ctrl_1', 'Ctrl_2', 'CD47_2', 'CD24_1']
L = []
for sample in samples:
    path_ch_matrix = os.path.join(path_data, sample)
    afm_raw = sc.read(os.path.join(path_ch_matrix, 'raw_RedeemV_untrimmed.h5ad'))
    L.append(afm_raw.obs[['Cas9']].assign(sample=sample))


##


# Supp Fig 20a: n cells, clones and SH
df = pd.concat(L)

# Fig 2a
fig, axs = plt.subplots(1,3,figsize=(4,1.8), sharey=True)

kwargs = {'orient':'h'}

df_ = df.groupby('sample')['Cas9'].nunique().to_frame().reset_index()
plu.bar(df_, x='Cas9', y='sample', ax=axs[0], color='#105D62', with_label=True, x_order=samples, kwargs=kwargs)
plu.format_ax(xlabel='n Cas9 clones', ylabel='', ax=axs[0], reduced_spines=True)

df_ = df.groupby('sample').size().to_frame('n cells').reset_index()
plu.bar(df_, x='n cells', y='sample', ax=axs[1], color='#105D62', with_label=True, x_order=samples, kwargs=kwargs)
plu.format_ax(xlabel='n cells', ylabel='', ax=axs[1], reduced_spines=True)

df_ = (
    df.groupby('sample')['Cas9']
    .apply(lambda x: - np.sum(x.value_counts(normalize=True) * np.log10(x.value_counts(normalize=True))))
    .to_frame('SH')
    .reset_index()
)
plu.bar(df_, x='SH', y='sample', ax=axs[2], color='#105D62', with_label=True, fmt="%.2f", x_order=samples, kwargs=kwargs)
plu.format_ax(xlabel='Shannon Entropy', ylabel='', ax=axs[2], reduced_spines=True)

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'Fig_20a.pdf'))


##