"""
Top MT-SNVs sets and associated plots.
"""

import os
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.dimred import *
from mito_utils.plotting_base import *
from mito_utils.diagnostic_plots import *
from mito_utils.embeddings_plots import *
from mito_utils.phylo_plots import *
matplotlib.use('macOSX')


##


# Set paths
sample ='MDA_PT'
job = 'tuning846006294d'
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro' 
path_results = os.path.join(path_main, f'results/MI_TO_bench/phylo_inference/final_trees/{sample}/{job}')
path_afm_raw = os.path.join(path_main, f'data/MI_TO_bench/AFMs/maegatk/{sample}/afm.h5ad')
path_coverage = os.path.join(path_main, f'data/MI_TO_bench/AFMs/maegatk/{sample}/coverage.txt.gz')
wobble_table = os.path.join(path_main, 'data/MI_TO_bench/miscellanea/formatted_table_wobble.csv')
weng_2024_ref = os.path.join(path_main, 'data/MI_TO_bench/miscellanea/weng2024_mut_spectrum_ref.csv')


##


# Read afm, raw: baseline preprocessing
afm_raw = sc.read(path_afm_raw)
afm_raw = filter_cells(afm_raw, cell_filter='filter2')
annotate_vars(afm_raw)
afm_raw = filter_baseline(afm_raw)

# Read processed afm
afm = sc.read(os.path.join(path_results, 'afm.h5ad'))
afm.uns['char_filter']

# UMAP, precomputed distances
reduce_dimensions(afm, metric='jaccard', bin_method='vanilla')


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


fig = plt.figure(figsize=(15,4.5))

ax = fig.add_subplot(1,4,1)
vars_AF_dist(afm_raw, ax=ax, color='#303030', alpha=.7, linewidth=.2)
vars_AF_dist(afm_raw[:,afm.var_names], ax=ax, color='#05A8B3', linewidth=.5, alpha=1)

ax = fig.add_subplot(1,4,2)
xticks = [1,2,4,10,30,90,300,1100]
plot_ncells_nAD(afm_raw, ax=ax,  xticks=xticks, c='#303030', s=2, alpha=.3)
plot_ncells_nAD(afm, ax=ax, c='#05A8B3', xticks=xticks, s=5, alpha=1, markeredgecolor='k')
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
fig.savefig(os.path.join(path_results, 'MT_SNVs.png'), dpi=500)
# plt.show()



##


# Mut profile
ref_df = pd.read_csv(weng_2024_ref, index_col=0)
fig = mut_profile(afm_raw.var_names, 
                  ref_df=ref_df, figsize=(5,2.5))
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'mut_profile_raw.png'), dpi=500)


##


# Colors
cmaps = {
    'GBC' : create_palette(afm.obs, 'GBC', sc.pl.palettes.default_102)
}

# Embeddings
gbcs = ['CTGCCTCCAACTAGCACC', 'TTACACGACCCGGCACGC', 'CATTGGATGGCGGTCACC']


afm.obs['GBC'].value_counts()

fig, axs = plt.subplots(2,3,figsize=(6,4.2))

for i,x in enumerate(gbcs):
    muts = []
    mut = compute_lineage_biases(afm, 'GBC', x).query('FDR<.1').sort_values('perc_in_target_lineage', ascending=False).index[0]
    afm.obs[mut] = afm[:,mut].X.A.flatten()
    s = f'GBC=="{x}"'
    draw_embeddings(afm.obs, cat='GBC', query=s, s=10, ax=axs[0,i],legend_kwargs={'colors':cmaps['GBC']}, 
                    axes_kwargs={'legend':False}, title=f'Clone {x[:5]}')
    axs[0,i].axis('off')  
    draw_embeddings(afm.obs.loc[lambda x: x[mut]<.001], cont=mut,  ax=axs[1,i], s=2,
                    cbar_kwargs={'palette':'mako_r', 'vmin':0, 'vmax':.01, 'label':'AF'}, axes_kwargs={'cbar':False})
    draw_embeddings(afm.obs.loc[lambda x: x[mut]>.001], cont=mut,  ax=axs[1,i], s=10, 
                    cbar_kwargs={'palette':'mako_r', 'vmin':0, 'vmax':.01, 'label':'AF'}, axes_kwargs={'cbar':False})
    axs[1,i].axis('off')

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'mut_embeddings.png'), dpi=500)


##


# All
fig, ax = plt.subplots(figsize=(4,4))
draw_embeddings(afm.obs, cat='GBC', s=50, ax=ax, legend_kwargs={'colors':cmaps['GBC']}, axes_kwargs={'legend':False})
format_ax(ax=ax, title=f'n cells: {afm.shape[0]}, n GBC: {afm.obs["GBC"].unique().size}', title_size=12)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'embeddings.png'), dpi=500)


##


# Trees
with open(os.path.join(path_results, 'annotated_tree.pickle'), 'rb') as f:
    tree = pickle.load(f)
tree_metrics = pd.read_csv(os.path.join(path_results, 'tree_metrics.csv'), index_col=0).set_index('metric')
_, mut_nodes, mutation_order = cut_and_annotate_tree(tree)
cmaps['MT_clone'] = create_palette(tree.cell_meta, 'MT_clone', sc.pl.palettes.godsnot_102)


##


fig, axs = plt.subplots(1,2,figsize=(8,4.5))

ax = axs[0]
plot_tree(tree, ax=ax, feature_internal_nodes='support', internal_node_subset=mut_nodes,
          cmap_internal_nodes=('Spectral_r',50,80), internal_node_kwargs={'markersize':5})
corr_distances = tree_metrics.loc["corr_distances", 'value']
supp_mut_nodes = np.median(get_supports(tree, subset=mut_nodes))
ax.set(title=f'Mut-assigned nodes support: {supp_mut_nodes}\nTree-MT-SNVs distances corr: {corr_distances:.2f}')
add_cbar(get_supports(tree), ax=ax, palette='Spectral_r', vmin=50, vmax=80, label='Support', layout='h3')

ax = axs[1]
plot_tree(tree, features=['GBC', 'MT_clone'], ax=ax, categorical_cmaps=cmaps, colorstrip_width=3)
ARI = tree_metrics.loc["ARI", 'value']
NMI = tree_metrics.loc["NMI", 'value']
ax.set(title=f'ARI: {ARI:.2f}, NMI: {NMI:.2f}')

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'phylo_main.png'), dpi=500)


##


fig, ax = plt.subplots(figsize=(4,4))

phylocorr = pd.read_csv(os.path.join(path_results, 'phylocorr.csv'), index_col=0)
ax.imshow(phylocorr.pivot_table(index='Var1', columns='Var2', values='Z'), vmin=-20,vmax=+20, cmap='inferno')
format_ax(ax=ax, xticks=[], yticks=[], xlabel='GBC', ylabel='GBC', xlabel_size=12, ylabel_size=12)
add_cbar(phylocorr['Z'], ax=ax, palette='inferno', vmin=-5, vmax=+5, label='Phylo-correlation', layout='outside')
fig.tight_layout()

fig.savefig(os.path.join(path_results, 'phylo_corr.png'), dpi=500)


##


# Main mut trees
fig, axs = plt.subplots(1,2,figsize=(15,8), gridspec_kw={'wspace': 0.3})

plot_tree(tree, ax=axs[0], 
    colorstrip_spacing=.000001, colorstrip_width=2,
    orient='down',
    features=['GBC', 'MT_clone']+mutation_order, layer='raw',
    categorical_cmaps=cmaps,
    feature_label_size=4, feature_label_offset=2
)
add_cbar(tree.layers['transformed'].values.flatten(), palette='mako', label='AF', ticks_size=8, label_size=9,
         vmin=.01, vmax=.1,
         ax=axs[0], layout=( (1.02,.3,.02,.2), 'right', 'vertical' ))

plot_tree(tree, ax=axs[1],
    colorstrip_spacing=.000001, colorstrip_width=2,
    orient='down',
    features=['GBC', 'MT_clone']+mutation_order, layer='transformed',
    feature_label_size=4, feature_label_offset=2,
    categorical_cmaps=cmaps,
    internal_node_subset=mut_nodes,
    internal_node_kwargs={'markersize':5, 'c':'darkred'}, show_internal=True
)
add_legend(label='Genotype', ax=axs[1], colors={'REF':'b', 'ALT':'r'}, loc='center left', bbox_to_anchor=(1,.4),
           ticks_size=8, artists_size=10, label_size=9)

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'mut_tree.png'))


##