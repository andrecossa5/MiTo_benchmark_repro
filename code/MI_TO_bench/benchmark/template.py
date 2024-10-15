"""
Coding game.
"""

import os
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.plotting_base import *
from mito_utils.dimred import *
from mito_utils.phylo import *
matplotlib.use('macOSX')

sample ='AML_clones'
path_afm = f'/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/data/MI_TO_bench/AFMs/freebayes/{sample}/afm.h5ad'
path_dbSNP = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/data/MI_TO_bench/miscellanea/dbSNP_MT.txt'
path_REDIdb = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/data/MI_TO_bench/miscellanea/REDIdb_MT.txt'
path_coverage = f'/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/data/MI_TO_bench/AFMs/freebayes/{sample}/coverage.txt.gz'
wobble_table = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/data/MI_TO_bench/miscellanea/formatted_table_wobble.csv'


afm = sc.read(path_afm)
afm = filter_cells(afm, cell_filter='None')
afm_raw = afm.copy()

afm, tree = filter_afm(
    afm,
    min_cell_number=10,
    lineage_column='GBC',
    filtering=None,
    # filtering_kwargs={
    #     'min_cov' : 10,
    #     'min_var_quality' : 30,
    #     'min_frac_negative' : .2,
    #     'min_n_positive' : 5,
    #     'af_confident_detection' : .01,
    #     'min_n_confidently_detected' : 2,
    #     'min_mean_AD_in_positives' : 1.5,       # 1.25,
    #     'min_mean_DP_in_positives' : 20
    # },
    binarization_kwargs={
        't_prob':.9, 't_vanilla':.0, 'min_AD':1, 'min_cell_prevalence':.2
    },
    bin_method='MI_TO',
    tree_kwargs={'metric':'custom_MI_TO_jaccard', 'solver':'NJ', 'ncores' : 8},
    path_dbSNP=path_dbSNP, 
    path_REDIdb=path_REDIdb,
    spatial_metrics=True,
    compute_enrichment=True,
    max_AD_counts=2,
    return_tree=True
)

# tree, _, _ = cut_and_annotate_tree(tree)
# stats = { k:v for k,v in afm.uns.items() }
# stats['n_MT_clone'] = tree.cell_meta['MT_clone'].nunique()
# stats['corr_dist'] = calculate_corr_distances(tree)

# reduce_dimensions(afm)


##


# Format coverage
cov = pd.read_csv(path_coverage, header=None)
cov.columns = ['pos', 'cell', 'n'] 
cov['cell'] = cov['cell'].map(lambda x: f'{x}_{sample}')
cov = cov.query('cell in @afm.obs_names')
cov['cell'] = pd.Categorical(cov['cell'], categories=afm.obs_names)
cov['pos'] = pd.Categorical(cov['pos'], categories=range(1,16569+1))
cov = cov.pivot_table(index='cell', columns='pos', values='n', fill_value=0)


##


from mito_utils.diagnostic_plots import *

fig = plt.figure(figsize=(15,4.5))

ax = fig.add_subplot(1,4,1)
vars_AF_dist(afm_raw, ax=ax, color='#303030', alpha=.7, linewidth=.2)
vars_AF_dist(afm_raw[:,afm.var_names], ax=ax, color='#05A8B3', linewidth=.5, alpha=1)

ax = fig.add_subplot(1,4,2)
xticks = [1,2,4,10,30,90,300,1100]
plot_ncells_nAD(afm_raw, ax=ax,  xticks=xticks, c='#303030', s=2, alpha=.3)
plot_ncells_nAD(afm, ax=ax, c='#05A8B3', xticks=xticks, s=3.5, alpha=.5)
format_ax(ax=ax, ylabel='Mean nAD / +cells', xlabel='n +cells', reduced_spines=True)

ax = fig.add_subplot(1,4,3, polar=True)
MT_coverage_polar(cov, var_subset=afm.var_names, ax=ax, kwargs_subset={'markersize':8, 'c':'#05A8B3'}, kwargs_main={'c':'#303030', 'linewidth':1.5, 'alpha':.7})

ax = fig.add_subplot(1,4,4)
ref_df = pd.read_csv(wobble_table, index_col=0)
ref_df['mut'] = ref_df['Position'].astype(str) + '_' + ref_df['Reference'] + '>' + ref_df['Variant']
df_ = ref_df.query('mut in @afm.var_names')['Symbol'].value_counts().to_frame('n')
bar(df_, 'n', ax=ax, c='#C0C0C0', edgecolor='k', s=.8)
format_ax(ax=ax, xticks=df_.index, rotx=90, ylabel='n MT-SNVs', xlabel='Gene', reduced_spines=True)

fig.subplots_adjust(bottom=.25, top=.8, left=.1, right=.9, wspace=.4)
fig.savefig('/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/results/MI_TO_bench/aa.png', dpi=500)

plt.show()


##


# Mut profile
ref_df = pd.read_csv('/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/data/MI_TO_bench/miscellanea/weng2024_mut_spectrum_ref.csv', index_col=0)
fig = mut_profile(afm.var_names, ref_df=ref_df, figsize=(5,2.5))
fig.tight_layout()
fig.savefig('/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/results/MI_TO_bench/bb.png', dpi=500)

fig.tight_layout()
plt.show()


##



# Embeddings
from mito_utils.embeddings_plots import *

gbcs = ['GTCGCTGTCCTGCTCCCG', 'CCCTTGCTTCCACTGTCC', 'CTGCCTCCAACTAGCACC']
reduce_dimensions(afm)

fig, axs = plt.subplots(2,3,figsize=(6,4.2))

for i,x in enumerate(gbcs):
    muts = []
    mut = compute_lineage_biases(afm, 'GBC', x).index[0]
    print(mut)
    afm.obs[mut] = afm[:,mut].X.A.flatten()
    s = f'GBC=="{x}"'
    draw_embeddings(afm.obs, cat='GBC', query=s, s=5, ax=axs[0,i], axes_kwargs={'legend':False}, title=f'Clone {x[:5]}')
    axs[0,i].axis('off')  
    draw_embeddings(afm.obs.loc[lambda x: x[mut]<.001], cont=mut,  ax=axs[1,i], s=2, cbar_kwargs={'palette':'mako_r', 'vmin':0, 'vmax':.05, 'label':'AF'}, axes_kwargs={'cbar':False})
    draw_embeddings(afm.obs.loc[lambda x: x[mut]>.01], cont=mut,  ax=axs[1,i], s=5, cbar_kwargs={'palette':'mako_r', 'vmin':0, 'vmax':.05, 'label':'AF'}, axes_kwargs={'cbar':False})
    axs[1,i].axis('off')

fig.tight_layout()
fig.savefig('/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/results/MI_TO_bench/cc.png', dpi=500)


##
