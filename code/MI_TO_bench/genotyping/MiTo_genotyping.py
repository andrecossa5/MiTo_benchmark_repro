"""
MiTo.
"""

import os
from mito_utils.utils import *
from mito_utils.preprocessing import *
from plotting_utils._plotting_base import *
matplotlib.use('macOSX')


##


# Args
path_afm = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/data/MI_TO_bench/AFMs/maegatk/MDA_PT/afm.h5ad'
path_dbSNP = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/data/MI_TO_bench/miscellanea/dbSNP_MT.txt'
path_REDIdb = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/data/MI_TO_bench/miscellanea/REDIdb_MT.txt'

afm = sc.read(path_afm)
afm = filter_cells(afm, cell_filter='filter_2')
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
        'min_mean_AD_in_positives' : 1.25,       # 1.25,
        'min_mean_DP_in_positives' : 5
    },
    binarization_kwargs={
        't_prob':.7, 't_vanilla':.0, 'min_AD':2, 'min_cell_prevalence':.05
    },
    bin_method='MiTo',
    path_dbSNP=path_dbSNP, 
    path_REDIdb=path_REDIdb,
    max_AD_counts=2
)

compute_distances(afm, precomputed=True)
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

# mut = '15876_T>C'
# mut = '9498_T>C'
mut = '9528_C>T'
ad = afm[:,mut].layers['AD'].A.flatten()
dp = afm[:,mut].layers['site_coverage'].A.flatten()
ad0, ad1 = get_components(ad[dp>0], dp[dp>0])


##


fig, ax = plt.subplots(figsize=(5,4))
colors = {'Observed':'k', 'C0':'#2EA867', 'C1':'#1D6C42'}
sns.kdeplot(ad[dp>0], ax=ax, color=colors['Observed'], bw_adjust=3, fill=True, alpha=.2)
sns.kdeplot(ad0, ax=ax, bw_adjust=5, color=colors['C0'], fill=True, alpha=.2)
sns.kdeplot(ad1, ax=ax, color=colors['C1'], bw_adjust=3, fill=True, alpha=.2)
add_legend(ax=ax, colors=colors, loc='upper right', label='Distribution', 
           bbox_to_anchor=(1,1), ticks_size=9, artists_size=10, label_size=10)
ax.set(xlabel='n of alternative UMIs')
ax.spines[['right', 'top']].set_visible(False)
fig.tight_layout()
plt.show()

#

df = genotype_mix(ad, dp, t_prob=.7, min_AD=2, debug=True)
df[df['geno_vanilla'] == df['geno_prob']]


##

