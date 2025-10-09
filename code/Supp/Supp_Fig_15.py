"""
Supp Fig 15
Re-analysis of Ludwig et. al, 2019 dataset: 3 GT-clones C9, D6, G10 in Supp Fig 2F. 
"""

import os
import numpy as np
import pandas as pd
import mito as mt
import matplotlib
import matplotlib.pyplot as plt
import plotting_utils as plu
from mito.tl.annotate import _get_muts_order
from scipy.io import mmread
from scipy.sparse import csr_matrix
from anndata import AnnData
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'data', 'ludwig_2019')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')

# Set visualization params
plu.set_rcParams({'figure.dpi':100})


##


# 1. Read data and format AFM --------------------------------------- # 

AD = pd.DataFrame(mmread(os.path.join(path_data, 'cellSNP.tag.AD.mtx')).toarray().T)
DP = pd.DataFrame(mmread(os.path.join(path_data, 'cellSNP.tag.DP.mtx')).toarray().T)

# Format AFM
AF = csr_matrix(np.divide(AD,(DP+.00000001)).values.astype(np.float64))
AD = csr_matrix(AD.values).astype(np.int64)
DP = csr_matrix(DP.values).astype(np.int64)
afm_original = AnnData(X=AF, 
              var=pd.DataFrame(index=[ f'MT-SNV {i}' for i in range(AD.shape[1]) ]),
              layers={'AD':AD, 'DP':DP},
              uns={'pp_method':'custom','scLT_system':'Smart-seq2'})

# 2. Explore data with MiTo defaults (tuned on MAESTER data) --------------------------------------- # 

# Try simple vanilla first
mt.pp.call_genotypes(afm_original)
mt.pp.compute_distances(afm_original, precomputed=True)

fig, ax = plt.subplots(figsize=(3,3))
mt.pl.heatmap_variants(afm_original, ax=ax, cmap='afmhot_r')
fig.tight_layout()
plt.show()             # Lot of noise and spurious counts, no clear structure, very high AF

# Annotate vars 
mt.pp.annotate_vars(afm_original)
afm_original.var

"""
Only 3 MT-SNVs with high AF > .1 in a substantial number of cells. 
Very high AF in +cells (all fixed variants, due to stringent filtering in original pubblication).
Differently from MAESTER data (UMI-based, consensus UMI calls in AD and DP counts), here MT-SNVs
are supported by high AD and DP raw read counts (no error correction by UMI collapsing).
Due to these features, we will tune MiTo parameters to be more stringent, and avoid 
spurious MT-SNVs calls.
"""

# Tuned genotyping
afm = afm_original.copy()
mt.pp.call_genotypes(afm, bin_method='MiTo', t_vanilla=.1, min_AD=10)

# Remove cells and variants with no calls (N.B., mimick MiTo filtering in mt.pp.filter_afm)
cells = afm.obs_names[(afm.layers['bin']==1).sum(axis=1).A1>0]
variants = afm.var_names[(afm.layers['bin']==1).sum(axis=0).A1>0]
afm = afm[cells,variants].copy()

# Compute distances
mt.pp.compute_distances(afm, precomputed=True)

# Given high noise and high AF of MT-SNVs in False Positive +cells, we
# re-adjust the parameters of clonal inference too. 
tree = mt.tl.build_tree(afm, precomputed=True)
model = mt.tl.MiToTreeAnnotator(tree)
model.clonal_inference(
    mut_enrichment_tresholds=[10],  # Higher MT-SNV enrichment treshold consider a node as clonal 
    af_treshold=.1                  # Higher pseudo-bulk AF treshold to call a variant in a clone
)


##


# 3. Supp Fig 15. Visualize AFM and phylogeny --------------------------------------- # 

fig, axs = plt.subplots(1,2,figsize=(7,4))

mt.pl.heatmap_variants(afm_original, ax=axs[0], cmap='afmhot_r', kwargs={'x_names_size':8})
plu.format_ax(ax=axs[0], rotx=90)
mt.pl.heatmap_variants(afm, tree=tree, ax=axs[1], cmap='afmhot_r', kwargs={'x_names_size':8})
plu.format_ax(ax=axs[1], rotx=90)

fig.subplots_adjust(wspace=0.5, top=.85, bottom=.22, left=.15, right=.85)
plu.save_best_pdf_quality(fig, figsize=(7,4), path=os.path.join(path_figures), name='Supp_Fig_15ab.pdf')


##


fig, ax = plt.subplots(figsize=(4,4))
mt.pl.plot_tree(
    tree, 
    ax=ax,
    colorstrip_width=2.3, 
    features=['MiTo clone'],
    internal_node_subset=tree.cell_meta['lca'].unique(), 
    show_internal=True,
    internal_node_kwargs={'markersize':5, 'c':'darkred'}
)
fig.tight_layout()
plu.save_best_pdf_quality(fig, figsize=(4,4), path=os.path.join(path_figures), name='Supp_Fig_15c.pdf')


##




