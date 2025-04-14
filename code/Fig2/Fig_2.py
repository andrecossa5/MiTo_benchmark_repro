"""
Benchmark output and MDA_PT example output.
"""

import os
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import mito as mt
import plotting_utils as plu
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_figures = os.path.join(path_main, 'results', 'figures', 'Fig2')
path_results = os.path.join(path_main, 'results', 'others', 'Fig2')


##


# Set params
plu.set_rcParams({'figure.dpi':150})

# Load clones colors
path_colors = os.path.join(path_main, 'data', 'general')
with open(os.path.join(path_colors, 'clones_colors_sc.pickle'), 'rb') as f:
    clones_colors = pickle.load(f)


##


# 1. MiTo bench experiment output viz --------------------------------------- # 

# Clonal composition: only final QC ('filter2' on MT modality) cells, before MT-SNVs filtering. 
samples = ['MDA_clones', 'MDA_PT', 'MDA_lung']
path_afms = os.path.join(path_main, 'data', 'general', 'AFMs', 'maegatk')

# Perform cell filtering (filter2), retrieve cell-clone assignment
L = []
for sample in samples:
    afm = sc.read(os.path.join(path_afms, f'afm_{sample}.h5ad'))
    afm = mt.pp.filter_cells(afm, cell_filter='filter2')
    L.append(afm.obs[['sample', 'GBC']])
df = pd.concat(L)
df['GBC'] = df['GBC'].astype('str')
qc_cells = df.index


##


fig, axs = plt.subplots(1,3,figsize=(4,1.5), sharey=True)

kwargs = {'orient':'h'}

df_ = df.groupby('sample')['GBC'].nunique().to_frame('n clones').reset_index()
plu.bar(df_, x='n clones', y='sample', ax=axs[0], color='#105D62', with_label=True, x_order=samples, kwargs=kwargs)
plu.format_ax(xlabel='n clones', ylabel='', ax=axs[0], reduced_spines=True)

df_ = df.groupby('sample').size().to_frame('n cells').reset_index()
plu.bar(df_, x='n cells', y='sample', ax=axs[1], color='#105D62', with_label=True, x_order=samples, kwargs=kwargs)
plu.format_ax(xlabel='n cells', ylabel='', ax=axs[1], reduced_spines=True)

df_ = (
    df.groupby('sample')['GBC']
    .apply(lambda x: - np.sum(x.value_counts(normalize=True) * np.log10(x.value_counts(normalize=True))))
    .to_frame('SH')
    .reset_index()
)
plu.bar(df_, x='SH', y='sample', ax=axs[2], color='#105D62', with_label=True, fmt="%.2f", x_order=samples, kwargs=kwargs)
plu.format_ax(xlabel='Shannon Entropy', ylabel='', ax=axs[2], reduced_spines=True)

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'clonal_numbers.pdf'))


##


fig, ax = plt.subplots(figsize=(4,1.5))

plu.bb_plot(df, 'sample', 'GBC', categorical_cmap='tab10', ax=ax)
fig.tight_layout()

fig.savefig(os.path.join(path_figures, 'clonal_composition.pdf'))


##


fig, axs = plt.subplots(1,3,figsize=(6,2),subplot_kw={'projection': 'polar'})

path_coverage = os.path.join(path_main, 'data', 'general', 'coverage', 'maegatk')

for i,sample in enumerate(samples):
    cov = pd.read_csv(os.path.join(path_coverage, f'{sample}_coverage.txt.gz'), header=None)
    cov.columns = ['pos', 'cell', 'n']
    cov['cell'] = cov['cell'].map(lambda x: f'{x}_{sample}')
    cov = cov.query('cell in @qc_cells')
    mt.pl.MT_coverage_by_gene_polar(cov, sample, ax=axs[i])

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'MT_coverage.pdf'))


##


# Legend genes
fig, ax = plt.subplots(figsize=(6,2))

df_mt = (
    pd.DataFrame(
        mt.ut.MAESTER_genes_positions, 
        columns=['gene', 'start', 'end']
    )
    .set_index('gene')
    .sort_values('start')   
)
cmap = { k:v for k,v in zip(df_mt.index, sc.pl.palettes.default_102[:df_mt.shape[0]])}
plu.add_legend(cmap, label='MT-genes', ax=ax, loc='center', frameon=True,
               bbox_to_anchor=(.5,.5), ncols=8, artists_size=5, ticks_size=6, label_size=8)
ax.axis('off')

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'MT_genes.pdf'))


##


# 2. MDA_PT example showcase --------------------------------------- # 

# Params
plu.set_rcParams({'figure.dpi':500})
# plt.rcParams

# See also explore output. ./results/others/explore_final_jobs
path_afms = os.path.join(path_main, 'data', 'general', 'AFMs', 'maegatk')
path_coverage = os.path.join(path_main, 'data', 'general', 'coverage', 'maegatk')
path_data = os.path.join(path_main, 'data', 'lineage_inference', 'UPMGA')

# Choose MT-SNVs space
sample = 'MDA_PT'
job_id = 'a4a3ca88dd'
path_data = os.path.join(path_data, sample, job_id)

# Load: afm, tree, metrics, phylocorr
afm_unfiltered = sc.read(os.path.join(path_afms, f'afm_{sample}.h5ad'))
afm = sc.read(os.path.join(path_data, 'afm_filtered.h5ad'))
with open(os.path.join(path_data, 'annotated_tree.pickle'), 'rb') as f:
    tree = pickle.load(f)


##


# Viz

# Cell and variant filtering
fig, axs = plt.subplots(1,3,figsize=(6,2))

# Cell filter
afm_unfiltered.obs['cell_filter'] = np.where(afm_unfiltered.obs_names.isin(afm.obs_names), 'Filtered', 'Discarded')
cmap = {'Discarded':'grey', 'Filtered':'#105D62'}

plu.scatter(afm_unfiltered.obs.query('cell_filter=="Discarded"'), 
            x='median_target_site_coverage', y='frac_target_site_covered', color=cmap['Discarded'],
            size=2, alpha=.3, ax=axs[0], kwargs={'zorder':1})
plu.scatter(afm_unfiltered.obs.query('cell_filter=="Filtered"'), 
            x='median_target_site_coverage', y='frac_target_site_covered', color=cmap['Filtered'],
            size=3, alpha=.3, ax=axs[0], kwargs={'zorder':2})

plu.format_ax(ax=axs[0], xlabel='Site coverage', ylabel='Fraction covered')
plu.add_legend(cmap, ax=axs[0], loc='lower right', bbox_to_anchor=(1,0))

# Var, filtered, AF spectrum
mt.pl.vars_AF_spectrum(afm_unfiltered[afm.obs_names], ax=axs[1], color='grey', linewidth=.1)
mt.pl.vars_AF_spectrum(afm_unfiltered[afm.obs_names,afm.var_names], ax=axs[1], color='#105D62', linewidth=1)

# Var, filtered, nAD, +cells
xticks = [1,2,5,15,50,200,800]
mt.pl.plot_ncells_nAD(afm_unfiltered[afm.obs_names], ax=axs[2], xticks=xticks, s=1, color='grey', alpha=.5)
mt.pl.plot_ncells_nAD(afm_unfiltered[afm.obs_names,afm.var_names], ax=axs[2], xticks=xticks, s=3, color='#105D62')

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'cell_variant_filters.pdf'))


##


# Mut signature
fig = mt.pl.mut_profile(afm.var_names, figsize=(4,1.8), legend_kwargs={'artists_size':6, 'ticks_size':6, 'label_size':6})
fig.savefig(os.path.join(path_figures, 'mut_profile.pdf'))


##


# Mut spatial organization
path_coverage = os.path.join(path_main, 'data', 'general', 'coverage', 'maegatk')
cov = pd.read_csv(os.path.join(path_coverage, f'{sample}_coverage.txt.gz'), header=None)
cov.columns = ['pos', 'cell', 'n']
cov['cell'] = cov['cell'].map(lambda x: f'{x}_{sample}')
cov = cov.query('cell in @afm.obs_names')

fig, ax = plt.subplots(figsize=(2,2),subplot_kw={'projection': 'polar'})
mt.pl.MT_coverage_polar(cov, var_subset=afm.var_names, ax=ax, 
                        xlabel_size=6, ylabel_size=8,
                        kwargs_subset={'markersize':2, 'c':'#105D62', 'marker':'o', 'alpha':.7})

ax.grid(True, linewidth=.3)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'mut_across_MT_genome.pdf'))


##


# Variants and distances heatmap
fig, axs = plt.subplots(1,2,figsize=(6,3))

mt.pl.heatmap_variants(afm, tree=tree, ax=axs[0], cmap='afmhot_r', kwargs={'x_names_size':3})
mt.pl.heatmap_distances(afm, tree=tree, ax=axs[1])

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'heatmaps.pdf'))


##


# Tree
fig, ax = plt.subplots(figsize=(4.5,4.5))
cmaps = { 
    'GBC' : clones_colors,
    'MiTo clone' : plu.create_palette(tree.cell_meta, 'MiTo clone', col_list=sc.pl.palettes.default_102)
}
clonal_nodes = tree.cell_meta['lca'].unique()[1:]
mt.pl.plot_tree(
    tree, 
    features=['GBC', 'MiTo clone'], 
    categorical_cmaps=cmaps, 
    ax=ax, 
    colorstrip_width=5,
    internal_node_subset=clonal_nodes,
    feature_internal_nodes='support',
    internal_node_kwargs={'markersize':5}
)
tree_stats = mt.tl.get_internal_node_stats(tree)
plu.add_cbar(
    tree_stats['support'], palette='Spectral_r', ax=ax, 
    vmin=.5, vmax=.8, label='Support', 
    layout=( (1-.27,.05,.22,.02), 'bottom', 'horizontal' )
)

fig.subplots_adjust(top=.9, bottom=.1, left=.1, right=.9)
fig.savefig(os.path.join(path_figures, f'tree.pdf'))
plt.close()


##


# Enriched variants

# UMAP
mt.pp.reduce_dimensions(afm, metric='jaccard')

# Main UMAP
fig, ax = plt.subplots(figsize=(4,4))
mt.pl.draw_embedding(afm, feature='GBC', ax=ax, categorical_cmap=clones_colors)
fig.tight_layout()
fig.savefig(os.path.join(path_figures, f'embeddings.pdf'))


##


# Select clones and variants
mapping = {
    'CACGATCCGTCGGCCCAG' : '15132_T>C',
    'TTACACGACCCGGCACGC' : '15660_C>T',
    'CTTGCCACTGATATCGAC' : '11814_T>C',
}

fig, axs = plt.subplots(3,2,figsize=(4.5,6))

for i,clone in enumerate(mapping):

    mut = mapping[clone]
    afm.obs['cat'] = np.where(afm.obs['GBC']==clone, clone, 'unassigned')
    _cmap = { clone : clones_colors[clone], 'unassigned' : 'white' }
    mt.pl.draw_embedding(afm, feature='cat', categorical_cmap=_cmap, ax=axs[i,0], outline=True)
    mt.pl.draw_embedding(afm, feature=mut, ax=axs[i,1], continuous_cmap='afmhot_r', 
                         legend=False, outline=True, kwargs={'colorbar_loc':None})
    axs[i,0].set(title=f'clone-{i+1}')
    axs[i,1].set(title=mut)
    plu.add_cbar(afm[:,mut].X.A.flatten(), palette='afmhot_r', ax=axs[i,1], label='AF')

fig.tight_layout()
fig.savefig(os.path.join(path_figures, f'enriched_muts.pdf'))


##

