"""
Supp Fig 20
- Re-analysis of Weng et. al, 2024 dataset: mouse batch 1.
- Pre-processing of 
  1: raw MT-SNVs calls from RedeemV
  2: filtered MT-SNVs calls after RedeemR v4 filtering 
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


# Set visualization params
plu.set_rcParams({'figure.dpi':120})


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'data', 'redeem_mouse')
path_results = os.path.join(path_main, 'results', 'others', 'Supp20')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')


# Pre-process MiTo AFM (add cas9 clones as columns), and gather cells meta
samples = ['DKO_1', 'CD47_1', 'Ctrl_1', 'Ctrl_2', 'CD47_2', 'CD24_1']

L = []
for sample in samples:

    path_ch_matrix = os.path.join(path_data, sample)
    
    # Cas9
    cas9 = sc.read(os.path.join(path_ch_matrix, 'cas9.h5ad')) 

    # Raw, untrimmed
    afm_raw = mt.io.make_afm(
        path_ch_matrix=path_ch_matrix, 
        scLT_system='RedeeM',
        pp_method='RedeemV',
        kwargs={'edge_trim': 0, 'treshold': 'Total'}
    )
    afm_raw.obs['Cas9'] = cas9.obs['cas9_clones']
    afm_raw.write(os.path.join(path_ch_matrix, 'raw_RedeemV_untrimmed.h5ad'))

    # Raw, Trimmed
    afm_raw = mt.io.make_afm(
        path_ch_matrix=path_ch_matrix, 
        scLT_system='RedeeM',
        pp_method='RedeemV',
        kwargs={'edge_trim': 5, 'treshold': 'Total'}
    )
    afm_raw.obs['Cas9'] = cas9.obs['cas9_clones']
    afm_raw.write(os.path.join(path_ch_matrix, 'raw_RedeemV_trimmed.h5ad'))

    ##

    # RedeemR, untrimmed
    afm_filtered = mt.io.make_afm(
        path_ch_matrix=path_ch_matrix, 
        scLT_system='RedeeM',
        pp_method='RedeemR',
        kwargs={'edge_trim': 0, 'treshold': 'Total'}
    )
    afm_filtered.obs['Cas9'] = cas9.obs['cas9_clones']
    afm_filtered.write(os.path.join(path_ch_matrix, 'filtered_RedeemR_untrimmed.h5ad'))
    # RedeemR, trimmed
    afm_filtered = mt.io.make_afm(
        path_ch_matrix=path_ch_matrix, 
        scLT_system='RedeeM',
        pp_method='RedeemR',
        kwargs={'edge_trim': 5, 'treshold': 'Total'}
    )
    afm_filtered.obs['Cas9'] = cas9.obs['cas9_clones']
    afm_filtered.write(os.path.join(path_ch_matrix, 'filtered_RedeemR_trimmed.h5ad'))


##

