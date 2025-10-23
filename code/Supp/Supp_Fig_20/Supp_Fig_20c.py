"""
Supp Fig 20c
- Re-analysis of Weng et. al, 2024 dataset: mouse batch 1.
- Visualize AOC, RF distances and internal node supports between in MT and Cas9 trees.
"""

import os
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import mito as mt
import matplotlib.pyplot as plt
import matplotlib
import plotting_utils as plu
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'results', 'others', 'Supp20')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')


##


# Extract sample and job_id couples
combos = []
for root, _, files in os.walk(os.path.join(path_data, 'trees_MT')):
    for file in files:
        if file.endswith('afm_filtered.h5ad'):
            sample = root.split('/')[-2]
            job_id = root.split('/')[-1]
            combos.append((sample, job_id))
            break

# AOCs_MT
AOCs_MT = []
AOCs_CAS = []
for sample, job_id in combos:

    # Read AFMs
    mito = sc.read(os.path.join(path_data, 'trees_MT', sample, job_id, 'afm_filtered.h5ad'))
    cas9 = sc.read(os.path.join(path_data, 'trees_Cas9', sample, job_id, 'afm_filtered.h5ad'))
    # Test same cells
    cells = list(set(mito.obs_names) & set(cas9.obs_names))
    assert len(cells) == len(set(mito.obs_names))
    # Re-order
    mito = mito[cells,:].copy()
    cas9 = cas9[cells,:].copy()
    # Extract distance matrices
    MT = mito.obsp['distances'].toarray()
    CAS = cas9.obsp['distances'].toarray()

    # Compute AOCs and related p-values
    aoc_mt, pvals_mt = mt.ut.AOC(MT, CAS, k=10, n_trials=1000)
    aoc_cas, pvals_cas = mt.ut.AOC(CAS, MT, k=10, n_trials=1000)
    df_mt = pd.DataFrame({'AOC':aoc_mt, 'pval':pvals_mt}, index=cells).assign(sample=sample, job_id=job_id)
    df_cas = pd.DataFrame({'AOC':aoc_cas, 'pval':pvals_cas}, index=cells).assign(sample=sample, job_id=job_id)
    AOCs_MT.append(df_mt)
    AOCs_CAS.append(df_cas)

# Concat AOCs_MT
df_aocs_mt = pd.concat(AOCs_MT).assign(direction='MT --> CAS9')
df_aocs_cas = pd.concat(AOCs_CAS).assign(direction='CAS9 --> MT')


##


# TBE supports
supports = []
for sample, job_id in combos:

    # Read trees
    with open(os.path.join(path_data, 'trees_MT', sample, job_id, 'annotated_tree.pickle'), 'rb') as f:
        tree_mt = pickle.load(f)
    with open(os.path.join(path_data, 'trees_Cas9', sample, job_id, 'annotated_tree.pickle'), 'rb') as f:
        tree_cas9 = pickle.load(f)

    supp_mt = mt.ut.get_internal_node_stats(tree_mt).assign(sample=sample, job_id=job_id, marker='MT')
    supp_cas9 = mt.ut.get_internal_node_stats(tree_cas9).assign(sample=sample, job_id=job_id, marker='Cas9')
    supports.append(supp_mt)
    supports.append(supp_cas9)

# Concat supports
df_supports = pd.concat(supports)


##


# RF distances
RFs = []
for sample, job_id in combos:

    # Read trees
    with open(os.path.join(path_data, 'trees_MT', sample, job_id, 'annotated_tree.pickle'), 'rb') as f:
        tree_mt = pickle.load(f)
    with open(os.path.join(path_data, 'trees_Cas9', sample, job_id, 'annotated_tree.pickle'), 'rb') as f:
        tree_cas9 = pickle.load(f)

    # RF calculations
    import cassiopeia as cs
    d, max_d = cs.critique.robinson_foulds(tree_cas9, tree_mt) 
    RFs.append(d / max_d)

# Concat RF distances
df_RF = (
    pd.DataFrame({'RF_distance':RFs}, index=[s for s, _ in combos])
    .sort_values('RF_distance')
    .reset_index(names='sample')
    .assign(mock='')
)


##


# Viz metrics distributions
fig, axs = plt.subplots(1,4,figsize=(4,1.8), width_ratios=[2,1.5,1.5,1.5])

# Supports
colors = plu.create_palette(df_supports, 'marker', palette=plu.darjeeling)
plu.dist(df_supports, 'support', by='marker', categorical_cmap=colors, ax=axs[0])
plu.add_legend(colors, ax=axs[0], loc='upper left', bbox_to_anchor=(-0.1,1.1), ticks_size=6, artists_size=5)
plu.format_ax(ax=axs[0], xlabel='TBE\nsupport', 
              title=f'Internal nodes',
              ylabel='Density', reduced_spines=True)

# AOC: MT --> Cas9
plu.dist(df_aocs_mt, 'AOC', color='#105D62', ax=axs[1])
plu.format_ax(ax=axs[1], xlabel='AOC\nMT>Cas9', 
              title=f'Cells',
              ylabel='', reduced_spines=True)

# AOC: Cas9 --> MT
plu.dist(df_aocs_cas, 'AOC', color='#105D62', ax=axs[2])
plu.format_ax(ax=axs[2], xlabel='AOC\nCas9>MT', 
              title=f'Cells',
              ylabel='', reduced_spines=True)

# RF distances
import seaborn as sns
colors = plu.create_palette(df_RF, 'sample', plu.darjeeling)
sns.stripplot(data=df_RF, x='mock', y='RF_distance', 
              palette=colors,
              hue='sample', dodge=False, ax=axs[3], legend=False)
# plu.dist(df_RF, 'RF_distance', color='#105D62', ax=axs[3])
plu.format_ax(ax=axs[3], xlabel='RF\ndistance', 
              title=f'Trees',
              ylabel='', reduced_spines=True)
axs[3].set_ylim((.965, 1.02))

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'Supp_Fig_20c.pdf'), dpi=300)
    

##