"""
Top MT-SNVs sets and associated plots.
"""

import os
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.dimred import *
from mito_utils.clustering import *
from mito_utils.phylo import *
from mito_utils.plotting_base import *
from mito_utils.diagnostic_plots import *
from mito_utils.embeddings_plots import *
from mito_utils.phylo_plots import *
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro' 
path_afms = os.path.join(path_main, 'data/MI_TO_bench/AFMs/maegatk')
path_results = os.path.join(path_main, 'results/MI_TO_bench/longitudinal')
path_meta = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/data/MDA/full_embs.csv'
path_dbSNP = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/data/MI_TO_bench/miscellanea/dbSNP_MT.txt'
path_REDIdb = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro/data/MI_TO_bench/miscellanea/REDIdb_MT.txt'


##


# Concat afms
PT = sc.read(os.path.join(path_afms, 'MDA_PT', 'afm.h5ad'))
lung = sc.read(os.path.join(path_afms, 'MDA_lung', 'afm.h5ad'))
afm = anndata.concat([PT, lung])
afm.uns['pp_method'] = 'maegatk'
afm.var['pos'] = afm.var_names.map(lambda x: x.split('_')[0]).astype(int)
afm.var['ref'] = afm.var_names.map(lambda x: x.split('_')[1].split('>')[0])
afm.var['alt'] = afm.var_names.map(lambda x: x.split('>')[1])
afm.obs['origin'] = np.where(afm.obs['sample']=='MDA_PT', 'PT', 'lung')

# Get cell states
meta = pd.read_csv(path_meta, index_col=0)
meta = meta.query('dataset=="NT_NT_2"')
meta.index = meta.index.map(lambda x: x.split('_')[0])
meta['new_cell_names'] = 'aa'
meta['new_cell_names'][meta['origin'] == 'PT'] = meta.index[meta['origin'] == 'PT'].map(lambda x: f'{x}_MDA_PT')
meta['new_cell_names'][meta['origin'] == 'lung'] = meta.index[meta['origin'] == 'lung'].map(lambda x: f'{x}_MDA_lung')
meta = meta.set_index('new_cell_names')
afm.obs = afm.obs.join(meta[['final_cell_state']])
afm.obs = afm.obs.rename(columns={'origin':'Origin', 'final_cell_state':'Cell state'})
afm.obs['Cell state'].loc[lambda x: x.isna()] = 'Undefined'
afm


##


# Filter cells
afm = filter_cells(afm, cell_filter='filter2')
afm = filter_afm(afm, 
    filtering='MI_TO', 
    filtering_kwargs={'af_confident_detection':.03, 'min_n_confidently_detected':2, 'min_mean_AD_in_positives':1.25},
    bin_method='MI_TO', 
    binarization_kwargs={'min_AD':2, 't_prob':.75, 'min_cell_prevalence':.1},
    max_AD_counts=2,
    lineage_column='GBC',
    min_cell_number=10,
    path_dbSNP=path_dbSNP, 
    path_REDIdb=path_REDIdb
)


##


# Cell states
tree.cell_meta.columns
tree.cell_meta.groupby('Origin')['GBC'].nunique()

##


# GBC 


##


# MT

##


# n MTs
fig, ax = plt.subplots(figsize=(4.5,4.5))
sns.kdeplot((afm.layers['AD'].A>0).sum(axis=1), 
            bw_adjust=1.5, ax=ax, fill=True, color='#00859A', linewidth=2.5, alpha=.5)
median = np.median((afm.layers['AD'].A>0).sum(axis=1))
ax.axvline(median, linestyle='--', c='k')
ax.set(xlabel='n MT-SNVs', title=f'Median n MT-SNVS per cell: {median}')
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'n MT-SNVs.png'), dpi=500)
plt.show()

# UMAP
reduce_dimensions(afm)

# All
fig, axs = plt.subplots(1,2,figsize=(8,4))
draw_embeddings(afm.obs, cat='Origin', s=10, ax=axs[0], legend_kwargs={'colors':{'PT':'#00AAC4', 'lung':'#C44F00'}}, axes_kwargs={'legend':False})
draw_embeddings(afm.obs, cat='GBC', s=10, ax=axs[1], axes_kwargs={'legend':False})
format_ax(ax=axs[1], title=f'n cells: {afm.shape[0]}, n GBC: {afm.obs["GBC"].unique().size}', title_size=12)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'PT_lung_embeddings.png'), dpi=500)
plt.show()



##


# Tree all
tree = build_tree(afm, metric='custom_MI_TO_jaccard', bin_method='MI_TO', solver='UPMGA', ncores=8,
                binarization_kwargs={'min_AD':2, 't_prob':.75, 'min_cell_prevalence':.1})
tree, mut_nodes, mutation_order = cut_and_annotate_tree(tree)


##


D = afm.obsp['distances'].A
order = leaves_list(linkage(D))
plt.imshow(1-D[np.ix_(order,order)], cmap='inferno')
plt.show()


from sklearn.metrics import normalized_mutual_info_score
custom_ARI(tree.cell_meta['MT_clone'], tree.cell_meta['GBC'])
normalized_mutual_info_score(tree.cell_meta['MT_clone'], tree.cell_meta['GBC'])
CI(tree).mean()
CI(tree).size


ci = CI(tree)

fig, ax = plt.subplots(figsize=(4,5))
order_gbc = tree.cell_meta['GBC'].value_counts().sort_values(ascending=False).index
order_mt = tree.cell_meta['MT_clone'].value_counts().sort_values(ascending=False).index
plot_heatmap(pd.crosstab(tree.cell_meta['MT_clone'], tree.cell_meta['GBC']).loc[order_mt,order_gbc], ax=ax, x_names_size=5, y_names_size=5)
fig.tight_layout()
plt.show()

# cells_ = tree.cell_meta['MT_clone'][tree.cell_meta['MT_clone']=='C12'].index
# afm.obs.loc[cells_,['median_target_site_coverage', 'median_untarget_site_coverage', 'frac_target_site_covered']].describe()
# afm.obs.loc[~tree.cell_meta.index.isin(cells_),['median_target_site_coverage', 'median_untarget_site_coverage', 'frac_target_site_covered']].describe()



##



# Colors
colors_gbc = create_palette(afm.obs, 'GBC', sc.pl.palettes.default_102)
colors_cell_states = create_palette(afm.obs, 'Cell state', sc.pl.palettes.vega_20)
colors_mt_clones = create_palette(tree.cell_meta, 'MT_clone', sc.pl.palettes.godsnot_102)
colors_mt_clones['Undefined'] = 'lightgrey'
colors_cell_states['Undefined'] = 'lightgrey'

# Tree plot
fig, ax = plt.subplots(figsize=(6,6))

# tree.collapse_mutationless_edges(True)
plot_tree(
    tree, ax=ax, 
    features=['Origin', 'GBC','MT_clone','Cell state']+mutation_order, 
    orient='down', layer='transformed',
    colorstrip_spacing=0.0000001,
    categorical_cmaps={
        'Origin':{'PT':'#00AAC4', 'lung':'#C44F00'}, 
        'GBC':colors_gbc,
        'MT_clone': colors_mt_clones,
        'Cell state':colors_cell_states,
    }, 
    feature_label_size=3,
    colorstrip_width=10, internal_node_subset=mut_nodes, internal_node_kwargs={'markersize':4, 'c':'r'}
)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'tree_all_prova.png'), dpi=500)
plt.show()

tree.character_matrix.describe()
CI(tree)


##


# Cell state legend
fig, ax = plt.subplots(figsize=(8,5))
add_legend(ax=ax, label='Cell state', colors=colors_cell_states, ncols=round(len(colors_cell_states)/4), bbox_to_anchor=(.5,.5), loc='center',
           ticks_size=8, label_size=10, artists_size=8)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'legend_cell_state.png'), dpi=500)


##


# Metastatic expansion clones.
early_clone = 'GTCGCTGTCCTGCTCCCG'
cells = afm.obs['GBC'].loc[lambda x: x==early_clone].index
early = afm[cells,:].copy()
# early = filter_afm(
#     early, 
#     filtering='MI_TO',
#     filtering_kwargs={'af_confident_detection':.1, 'min_n_confidently_detected':2, 'min_mean_AD_in_positives':1.25},
#     bin_method='vanilla',#'MI_TO',
#     binarization_kwargs={'min_AD':2, 't_prob':.5, 't_vanilla':.02, 'min_cell_prevalence':.2}
# )
# 
annotate_vars(early, overwrite=True)
early = early[:,early.var['n5']>5]
early = early[(early.X.A>.05).any(axis=1),:].copy()
compute_distances(early, metric='custom_MI_TO_jaccard', bin_method='MI_TO', 
                  binarization_kwargs={'t_vanilla':.95, 'min_cell_prevalence':.1, 'min_AD':2})
np.sum(early.layers['bin']>0)
early_tree = build_tree(early, precomputed=True, solver='UPMGA')

fig, ax = plt.subplots(figsize=(4,4))
order = leaves_list(linkage(early.obsp['distances'].A))
ax.imshow(early.obsp['distances'].A[np.ix_(order,order)], cmap='inferno_r')
plt.show()

fig, ax = plt.subplots(figsize=(4,4))
D = pairwise_distances(early.layers['bin'].A.T, metric='jaccard')
order = leaves_list(linkage(D))
ax.imshow(D[np.ix_(order,order)], cmap='inferno_r')
plt.show()

early_tree, mut_nodes, mut_order = cut_and_annotate_tree(early_tree, n_clones=3)
for node in early_tree.internal_nodes:
    early_tree.set_attribute(node, 'aa', 1)

fig, axs = plt.subplots(1,2,figsize=(8,4), gridspec_kw={'wspace':.3})
plot_tree(
    early_tree, ax=axs[0], features=['Origin']+mut_order, 
    orient='down', layer='raw',
    categorical_cmaps={
        'Origin':{'PT':'#00AAC4', 'lung':'#C44F00'}, 
    }, 
    colorstrip_width=2, internal_node_subset=mut_nodes, 
    feature_internal_nodes='aa', 
    internal_node_kwargs={'markersize':5, 'c':'r'}
)
plot_tree(
    early_tree, ax=axs[1], features=['Origin']+mut_order, 
    orient='down', layer='transformed',
    categorical_cmaps={
        'Origin':{'PT':'#00AAC4', 'lung':'#C44F00'}
    }, 
    colorstrip_width=2, internal_node_subset=mut_nodes, 
    feature_internal_nodes='aa', 
    internal_node_kwargs={'markersize':5, 'c':'r'}
)
fig.tight_layout()
# plt.show()
fig.savefig(os.path.join(path_results, 'tree_early.png'), dpi=500)


##


# Muts characterization
df_clones = (
    afm.obs.groupby(['Origin', 'GBC'])
    .size().to_frame('n').reset_index()
    .pivot(index='GBC', columns='Origin', values='n')
    .fillna(0)
    .query('PT>10 and lung>10')
    .assign(n_mean= lambda x: (x['PT']+x['lung']) / 2)
    .sort_values('n_mean', ascending=False)
)

clone_muts = {}
for clone in df_clones.index:

    clone_muts[clone] = {}

    enriched_pt_muts = (
        compute_lineage_biases(afm[afm.obs['Origin']=='PT'], 'GBC', clone, bin_method='vanilla', binarization_kwargs={'min_AD':2}, alpha=0.05)
        .query('perc_in_target_lineage>=.75 and FDR<=.1').index.to_list()
    )
    clone_muts[clone]['Enriched PT'] = len(enriched_pt_muts)
    enriched_lung_muts = (
        compute_lineage_biases(afm[afm.obs['Origin']=='lung'], 'GBC', clone, bin_method='vanilla', binarization_kwargs={'min_AD':2}, alpha=0.05)
        .query('perc_in_target_lineage>=.75 and FDR<=.1').index.to_list()
    )
    clone_muts[clone]['Enriched lung'] = len(enriched_lung_muts)
    clone_muts[clone]['Enriched both'] = len( set(enriched_pt_muts) & set(enriched_lung_muts) )

    pt_cells = afm.obs.query('GBC==@clone and Origin=="PT"').index
    PT = filter_afm(afm, filtering=None, cells=pt_cells, bin_method='vanilla', binarization_kwargs={'min_AD':2})
    test = np.sum(PT.layers['bin'].A==1, axis=0)>0
    pt_vars = PT.var_names[test]
    clone_muts[clone]['PT'] = np.sum(test)

    lung_cells = afm.obs.query('GBC==@clone and Origin=="lung"').index
    lung = filter_afm(afm, filtering=None, cells=lung_cells, bin_method='vanilla', binarization_kwargs={'min_AD':2})
    test = np.sum(lung.layers['bin'].A==1, axis=0)>0
    lung_vars = lung.var_names[test]
    clone_muts[clone]['Lung'] = np.sum(test)

    enriched_both_vars = set(enriched_pt_muts) & set(enriched_lung_muts)
    acquired_vars = set(lung_vars) - set(pt_vars)
    clone_muts[clone]['Acquired'] = len(acquired_vars)
    lost_vars = set(pt_vars) - set(lung_vars)
    clone_muts[clone]['Lost'] = len(lost_vars)
    
    clone_muts[clone]['_Enriched both'] = np.nanmean(
        np.where(afm[pt_cells,list(enriched_both_vars)].X.A>0, afm[pt_cells,list(enriched_both_vars)].X.A, np.nan)
    )
    clone_muts[clone]['_Acquired'] = np.nanmean(
        np.where(afm[pt_cells,list(acquired_vars)].X.A>0, afm[pt_cells,list(acquired_vars)].X.A, np.nan)
    )
    clone_muts[clone]['_Lost'] = np.nanmean(
        np.where(afm[pt_cells,list(lost_vars)].X.A>0, afm[pt_cells,list(lost_vars)].X.A, np.nan)
    )

df_clones = pd.DataFrame(clone_muts).T


## 


# Viz
fig, axs = plt.subplots(1,2,figsize=(8,4.5))

df_ = df_clones.loc[:,~df_clones.columns.str.startswith('_')].reset_index(names='GBC').melt(id_vars='GBC')
order = ['PT', 'Lung', 'Enriched PT', 'Enriched lung', 'Enriched both', 'Acquired', 'Lost']
strip(df_, x='variable', y='value', c='k', s=5, ax=axs[0], order=order)
format_ax(ax=axs[0], title='', rotx=45, ylabel='n MT-SNVs', reduced_spines=True)

df_ = df_clones.loc[:,df_clones.columns.str.startswith('_')].reset_index(names='GBC').melt(id_vars='GBC')
df_['variable'] = df_['variable'].map(lambda x: x[1:])
order = ['Enriched both', 'Acquired', 'Lost']
strip(df_, x='variable', y='value', c='k', s=5, ax=axs[1], order=order)
format_ax(ax=axs[1], title='', rotx=45, ylabel='Mean AF +cells', reduced_spines=True)

fig.suptitle('MT-SNVs after ~1 month')
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'MT_dynamics.png'), dpi=500)


##
