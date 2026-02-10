"""
Supp Fig 20d
- Re-analysis of Weng et. al, 2024 dataset: mouse batch 1.
- Visualize selected tree couples.
"""

import os
import pickle
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

# Set visualization params
plu.set_rcParams({'figure.dpi':120})


##


# Viz

# Extract sample and job_id couples
combos = []
for root, _, files in os.walk(os.path.join(path_data, 'trees_MT')):
    for file in files:
        if file.endswith('afm_filtered.h5ad'):
            sample = root.split('/')[-2]
            job_id = root.split('/')[-1]
            combos.append((sample, job_id))
            break

# Plot trees
fig, axs = plt.subplots(2,4,figsize=(12,6.5))

for i, (sample, job_id) in enumerate(combos):

    # Get axes
    ax_mito = axs[0,i]
    ax_cas9 = axs[1,i]

    # Read trees
    with open(os.path.join(path_data, 'trees_MT', sample, job_id, 'annotated_tree.pickle'), 'rb') as f:
        tree_mt = pickle.load(f)
    with open(os.path.join(path_data, 'trees_Cas9', sample, job_id, 'annotated_tree.pickle'), 'rb') as f:
        tree_cas9 = pickle.load(f)

    # Harmonize clone column names
    tree_cas9.cell_meta = (
        tree_cas9.cell_meta
        .rename(columns={'MiTo_clones':'MiTo clone', 'cas9_clones':'Cas9'})
    )

    # Plot
    mt.pl.plot_tree(tree_mt, ax=ax_mito, features=['Cas9', 'MiTo clone'], colorstrip_width=2.5)
    plu.format_ax(ax=ax_mito, title=f'{sample} - MT')
    mt.pl.plot_tree(tree_cas9, ax=ax_cas9, features=['Cas9', 'MiTo clone'], colorstrip_width=4)
    plu.format_ax(ax=ax_cas9, title=f'{sample} - Cas9')

fig.tight_layout()
plu.save_best_pdf_quality(
    fig, figsize=(12,6.5), 
    path=path_figures,
    name='Supp_Fig_20d.pdf'
)