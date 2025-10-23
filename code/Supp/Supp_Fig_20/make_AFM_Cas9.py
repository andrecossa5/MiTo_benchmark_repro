"""
Supp Fig 20
- Re-analysis of Weng et. al, 2024 dataset: mouse batch 1.
- Subset Cas9 AFMs using final filtered AFMs from RedeemR nf-MiTo INFER after TUNE. 
"""

import os
import pickle
import scanpy as sc

##

# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'data', 'redeem_mouse')
path_results = os.path.join(path_main, 'results', 'others', 'Supp20')

##

# Subset AFM and add MiTo clones to original, demultiplexed Cas9 AFMs

# Only for dataset in the end
samples = ['Ctrl_1', 'Ctrl_2', 'CD47_2', 'CD24_1'] 

for sample in samples:

    # Retrieve path to filtered MiTo AFM, and its tree
    for root, _, files in os.walk(os.path.join(path_results, 'trees_MT', sample)):
        if 'annotated_tree.pickle' in files:
            job_id = root.split('/')[-1]
            path_filtered_mt = os.path.join(root, 'annotated_tree.pickle')
            break

    # Read MT tree
    with open(path_filtered_mt, 'rb') as f:
        mt_tree = pickle.load(f)

    # Read original Cas9 AFM
    cas9_original = sc.read(os.path.join(path_data, sample, 'cas9.h5ad'))

    # Subset Cas9 AFM using cells in final filtered MT tree, add MiTo clones for plotting (later)
    cas9 = cas9_original[mt_tree.leaves,:].copy()
    cas9.obs['MiTo_clones'] = mt_tree.cell_meta['MiTo clone']

    # Write out subsetted Cas9 AFM
    cas9.write(os.path.join(path_results, 'cas9_filtered', f'{sample}_{job_id}_filtered.h5ad'))


##