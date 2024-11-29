"""
MiTo.
"""

import os
from itertools import product
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.plotting_base import *
from mito_utils.metrics import *
from mito_utils.dimred import *
from mito_utils.phylo import *
matplotlib.use('macOSX')


##


# Args
os.chdir('/Users/IEO5505/Desktop/MI_TO/mito_utils/scratch')
path_dbSNP = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/data/MI_TO_bench/miscellanea/dbSNP_MT.txt'
path_REDIdb = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/data/MI_TO_bench/miscellanea/REDIdb_MT.txt'

afm = sc.read('afm.h5ad')
afm = filter_cells(afm, cell_filter='None')
afm_raw = afm.copy()

afm = filter_afm(
    afm,
    min_cell_number=10,
    lineage_column='GBC',
    filtering='MiTo',
    filtering_kwargs={
        'min_cov' : 10,
        'min_var_quality' : 30,
        'min_frac_negative' : .2,
        'min_n_positive' : 5,
        'af_confident_detection' : .03,
        'min_n_confidently_detected' : 2,
        'min_mean_AD_in_positives' : 1.5,       # 1.25,
        'min_mean_DP_in_positives' : 5
    },
    binarization_kwargs={
        't_prob':.9, 't_vanilla':.0, 'min_AD':2, 'min_cell_prevalence':.05
    },
    bin_method='MiTo',
    path_dbSNP=path_dbSNP, 
    path_REDIdb=path_REDIdb,
    max_AD_counts=2
)

compute_distances(afm, precomputed=True)
tree = build_tree(afm, precomputed=True)
tree, _,_ = MiToTreeAnnotator(tree)
df = tree.cell_meta.dropna()
ari = custom_ARI(df['GBC'], df['MT_clone'])
nmi = normalized_mutual_info_score(df['GBC'], df['MT_clone'])
ari
nmi

calculate_corr_distances(tree)

np.sum(afm.layers['AD'].A, axis=0)
np.sum(afm.layers['site_coverage'].A>0, axis=0)


##


# Muts
var_order = np.argsort((afm.layers['AD'].A).mean(axis=0))
afm.var_names[var_order][-10:]

mut = '9441_C>T'
ad = afm[:,mut].layers['AD'].A.flatten()
dp = afm[:,mut].layers['site_coverage'].A.flatten()


# Dist fit
ad_th, _ = fit_mixbinom(ad[dp>0],dp[dp>0])

fig, axs = plt.subplots(1,5,figsize=(10,2), sharey=True)

for i,mut in enumerate(afm.var_names[var_order][-5:]):
    ax = axs.ravel()[i]
    ad = afm[:,mut].layers['AD'].A.flatten()
    dp = afm[:,mut].layers['site_coverage'].A.flatten()
    ad_th, _ = fit_mixbinom(ad[dp>0],dp[dp>0])
    sns.kdeplot(ad, ax=ax, color='k')
    sns.kdeplot(ad_th, ax=ax, color='r')
    
fig.tight_layout()
plt.show()


# Mixbinom components
var_order = np.argsort((afm.layers['AD'].A).mean(axis=0))
afm.var_names[var_order][-10:]

mut = '15876_T>C'
ad = afm[:,mut].layers['AD'].A.flatten()
dp = afm[:,mut].layers['site_coverage'].A.flatten()

ad0, ad1 = get_components(ad[dp>0], dp[dp>0])


##

afm.var

fig, ax = plt.subplots(figsize=(4,4))

colors = {'Observed':'k', 'C0':'orange', 'C1':'darkred'}
sns.kdeplot(ad[dp>0], ax=ax, bw_adjust=3, color=colors['Observed'], fill=True, alpha=.2)
sns.kdeplot(ad0, ax=ax, bw_adjust=3, color=colors['C0'], fill=True, alpha=.2)
sns.kdeplot(ad1, ax=ax, color=colors['C1'], fill=True, alpha=.2)
add_legend(ax=ax, colors=colors, loc='upper right', label='Distribution', bbox_to_anchor=(1,1), ticks_size=9, artists_size=10, label_size=10)
fig.tight_layout()
plt.show()

#

df = genotype_mix(ad, dp, t_prob=.7, min_AD=2, debug=True)


df[df['geno_vanilla'] != df['geno_prob']]


##

