"""
Discretization. NB and tresholds.
"""

import os
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.genotyping import *
from mito_utils.genotyping import _genotype_mix
from mito_utils.dimred import *
from mito_utils.plotting_base import *
from mito_utils.embeddings_plots import *
from mito_utils.diagnostic_plots import *
from mito_utils.phylo import *
from mito_utils.phylo_plots import *


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro'
path_data = os.path.join(path_main, 'data', 'MI_TO_bench')
path_results = os.path.join(path_main, 'results',  'MI_TO_bench')
make_folder(path_results, 'discretization', overwrite=False)

# Params
sample = 'MDA_clones'
filtering = 'MI_TO'

# Make an annotated AF matrix for the sample
path_afm = os.path.join(path_data, 'AFMs', f'{sample}.h5ad')
path_meta = os.path.join(path_data, 'cells_meta.csv')
afm = read_one_sample(path_afm, path_meta, sample=sample, nmads=5, mean_coverage=25)

# Read and filter
a, report = filter_cells_and_vars(
    afm, filtering=filtering, 
    max_AD_counts=3, af_confident_detection=.01, min_cell_number=10,
    lineage_column='GBC', compute_enrichment=True,
    path_dbSNP=os.path.join(path_data, 'miscellanea', 'dbSNP_MT.txt'),
    path_REDIdb=os.path.join(path_data, 'miscellanea', 'REDIdb_MT.txt'),
)

(a.var.loc[:,a.var.columns.str.startswith('FDR')]<=.1).sum(axis=0)

AD, DP, _ = get_AD_DP(a)
AD = AD.A.T
DP = DP.A.T

# VAF spectrum
fig, ax = plt.subplots(figsize=(4.5,4.5))
vars_AF_dist(a, ax=ax, color='darkred', linewidth=1, linestyle='--')
fig.tight_layout()
plt.show()

(a.X>0).sum(axis=1)                 # n vars (per cell)
(a.X>0).sum(axis=0) / a.shape[0]    # Prop +cells (per var)
a.X.mean(axis=0)                    # Mean AF vars
AD.mean(axis=0)                     # Mean AD vars
DP.mean(axis=0)                     # Mean DP vars

fig, ax = plt.subplots(figsize=(4.5,4.5))
ax.plot((a.X>0).sum(axis=0), np.ma.masked_less_equal(AD, 0).mean(axis=0), 'ko')
format_ax(ax=ax, xlabel='n +cells', ylabel='Mean n of UMIs for ALT in +cells')
fig.tight_layout()
plt.show()




##


# P-P plot
order = np.argsort((a.X>0).sum(axis=0))
i_ = int(np.round(order.size / 2)-1)
variants = a.var_names[order[::-1][:2]].to_list() + a.var_names[order[i_:(i_+2)]].to_list() + a.var_names[order[:2]].to_list()
dists = ['binom', 'nbinom', 'betabinom', 'mixbinom']

fig, axs = plt.subplots(nrows=4, ncols=len(variants), figsize=(15,11), sharex=True, sharey=True)

i = 0
for row, dist in enumerate(dists):
    for variant in variants:

        idx = np.where(a.var_names==variant)[0][0]
        ad = AD[:,idx]
        dp = DP[:,idx]
        wt = dp - ad

        if dist == 'binom':
            ad_th, _ = fit_binom(ad, dp)  
        elif dist == 'nbinom':
            ad_th, _ = fit_nbinom(ad, dp)  
        elif dist == 'betabinom':
            ad_th, _ = fit_betabinom(ad, dp)    
        elif dist == 'mixbinom':
            ad_th, _ = fit_mixbinom(ad, dp)          

        # CDFs
        X = np.concatenate((ad, ad_th))
        bins = np.linspace(np.min(X)-0.5, np.max(X)+0.5)
        ad_counts, _ = np.histogram(ad, bins=bins)
        ad_th_counts, _ = np.histogram(ad_th, bins=bins) 
        empirical_cdf = ad_counts.cumsum() / ad_counts.sum()
        theoretical_cdf = ad_th_counts.cumsum() / ad_th_counts.sum()

        # Stats variables
        mean_af = (ad/(dp+0.00001)).mean()
        freq =  ((ad/(dp+0.00001))>0).sum()

        ax = axs.ravel()[i]
        ax.plot(theoretical_cdf, empirical_cdf, 'o-')
        ax.plot([0,1], [0,1], 'r--')
        corr, p = stats.pearsonr(empirical_cdf, theoretical_cdf)
        if dist == 'binom':
            ax.set(title=f'{variant} \n Mean af: {mean_af:.2f}, n+: {int(freq)} \n Corr: {corr:.2f}; p {p:.2f}', ylabel='Empirical CDF', xlabel=f'{dist} CDF')
        else:
            ax.set(title=f'Corr: {corr:.2f}; p {p:.2f}', ylabel='Empirical CDF', xlabel=f'{dist} CDF')
        
        # Update figure count
        i += 1

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'AD_counts_dist_CDF.png'), dpi=400)


##


# BIC and AD fits
fig, axs = plt.subplots(nrows=1, ncols=len(variants), figsize=(19,4))

dists_ = ['Observed'] + dists
colors = {k:v for k,v in zip(dists_, ten_godisnot)}

for i,variant in enumerate(variants):

    ax = axs.ravel()[i]
    pos_y = .9
    delta = .08

    for row, dist in enumerate(dists_):

        idx = np.where(a.var_names==variant)[0][0]
        ad = AD[:,idx]
        dp = DP[:,idx]
        wt = dp - ad

        if dist == 'binom':
            x, d = fit_binom(ad, dp)  
        elif dist == 'nbinom':
            x, d = fit_nbinom(ad, dp)  
        elif dist == 'betabinom':
            x, d = fit_betabinom(ad, dp)    
        elif dist == 'mixbinom':
            x, d = fit_mixbinom(ad, dp)      
        else:
            x = ad   
            d = None 

        # Stats variables
        bic = d['BIC'] if isinstance(d, dict) else None
        L = d['L'] if isinstance(d, dict) else None
        mean_af =  (ad/(dp+0.00001)).mean()
        freq =  ((ad/(dp+0.00001))>0).sum()

        sns.kdeplot(x, ax=ax, color=colors[dist], fill=True, alpha=.3)
        ax.set(title=f'{variant} \n Mean af: {mean_af:.2f}, n+: {int(freq)}', ylabel='Density', xlabel=f'n ALT')
        if d is not None:
            ax.text(.28, pos_y, f'{dist} BIC: {bic:.2f}', transform=ax.transAxes, size=9)
            pos_y -= delta
    
    if i == 2:
        add_legend('Distribution', colors=colors, ax=ax, bbox_to_anchor=(1.15,-.3), loc='center', ticks_size=10, ncols=len(dist), label_size=11)

fig.subplots_adjust(left=.05, right=.95, top=.8, bottom=.25, wspace=.4)
fig.savefig(os.path.join(path_results, 'AD_counts_dists.png'), dpi=300)


##


# 2 components mixture of binomials
variant = variants[0]
colors = { 'observed': 'k', 'c0' : '#217043', 'c1' : '#bd4c0f' }

fig, axs = plt.subplots(nrows=1, ncols=len(variants), figsize=(19,3.7))

for i,variant in enumerate(variants):

    idx = np.where(a.var_names==variant)[0][0]
    ad = AD[:,idx]
    dp = DP[:,idx]
    wt = dp - ad

    np.random.seed(1234)
    model = MixtureBinomial(n_components=2, tor=1e-20)
    model.fit((ad, dp), max_iters=500, early_stop=True)

    # Access the estimated parameters
    ps = model.params[:2]
    pis = model.params[2:]
    idx1 = np.argmax(ps)
    idx0 = 0 if idx1 == 1 else 1
    p1 = ps[idx1]
    p0 = ps[idx0]
    pi1 = pis[idx1]
    pi0 = pis[idx0]
    d = {'p':[p0,p1], 'pi':[pi0,pi1]}

    ad0_th = dp * p0
    ad1_th = dp * p1

    ax = axs.ravel()[i]
    sns.kdeplot(ad, ax=ax, color="k", fill=True, alpha=.3)
    sns.kdeplot(ad0_th, ax=ax, color="#217043", fill=True, alpha=.3)
    sns.kdeplot(ad1_th, ax=ax, color="#bd4c0f", fill=True, alpha=.3)
    ax.set(title=f'{variant} \n Mean af: {mean_af:.2f}, n+: {int(freq)}', ylabel='Density', xlabel=f'n ALT')
    if i==0:
        add_legend('Distribution', colors=colors, ax=ax, bbox_to_anchor=(1,1), loc='upper right', ticks_size=10, label_size=10, artists_size=10)

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'AD_counts_dist_bibonm_mixture.png'), dpi=300)


##







