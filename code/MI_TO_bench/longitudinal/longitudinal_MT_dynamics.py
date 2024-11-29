"""
Top MT-SNVs sets and associated plots.
"""

import os
from itertools import chain
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.plotting_base import *
from mito_utils.diagnostic_plots import *
from mito_utils.phylo_plots import *
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro' 
path_afms = os.path.join(path_main, 'data/MI_TO_bench/AFMs/maegatk')
path_variants = os.path.join(path_main, 'results/MI_TO_bench/phylo_inference')
path_results = os.path.join(path_main, 'results/MI_TO_bench/longitudinal')
path_meta = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/data/MDA/full_embs.csv'


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


# Variants
vars_PT = set(pd.read_csv(os.path.join(path_variants, 'MDA_PT_top_10_var_subset.csv'), index_col=0)['var'].unique())
vars_lung = set(pd.read_csv(os.path.join(path_variants, 'MDA_lung_top_10_var_subset.csv'), index_col=0)['var'].unique())
variants = list(vars_lung | vars_PT)

# Filter AFM
bin_method = 'MiTo'
binarization_kwargs = {'min_AD':2, 't_prob':.75, 'min_cell_prevalence':.1}
afm = filter_cells(afm, cell_filter='filter2')
variants = afm.var_names[afm.var_names.isin(variants)]
afm = filter_afm(afm, filtering=None, variants=variants, bin_method=bin_method, binarization_kwargs=binarization_kwargs, min_cell_number=10, lineage_column='GBC')


##


# Muts characterization
df_clones = (
    afm.obs.groupby(['Origin', 'GBC'])
    .size().to_frame('n').reset_index()
    .pivot(index='GBC', columns='Origin', values='n')
    .fillna(0)
    .query('PT>=10 and lung>=10')
    .assign(n_mean= lambda x: (x['PT']+x['lung']) / 2)
    .sort_values('n_mean', ascending=False)
)

clone_muts = {}
enriched_muts = {}
selected_muts = {}

for clone in df_clones.index:

    clone_muts[clone] = {}
    enriched_pt_muts = (
        compute_lineage_biases(afm[afm.obs['Origin']=='PT'], 'GBC', clone)
        .query('prevalence>=.01 and perc_in_target_lineage>=0.75 and FDR<=.1').index.to_list()
    )
    clone_muts[clone]['Enriched PT'] = len(enriched_pt_muts)
    enriched_lung_muts = (
        compute_lineage_biases(afm[afm.obs['Origin']=='lung'], 'GBC', clone)
        .query('prevalence>=.01 and perc_in_target_lineage>=0.75 and FDR<=.1').index.to_list()
    )
    clone_muts[clone]['Enriched lung'] = len(enriched_lung_muts)
    clone_muts[clone]['Enriched both'] = len( set(enriched_pt_muts) & set(enriched_lung_muts) )
    enriched_muts[clone] = set(enriched_pt_muts) | set(enriched_lung_muts)

    pt_cells = afm.obs.query('GBC==@clone and Origin=="PT"').index
    PT = afm[pt_cells,:].copy()
    test = np.sum(PT.layers['bin'].A==1, axis=0)>=5
    pt_vars = PT.var_names[test]
    clone_muts[clone]['PT'] = np.sum(test)

    lung_cells = afm.obs.query('GBC==@clone and Origin=="lung"').index
    lung = afm[lung_cells,:].copy()
    test = np.sum(lung.layers['bin'].A==1, axis=0)>=5
    lung_vars = lung.var_names[test]
    clone_muts[clone]['Lung'] = np.sum(test)
    selected_muts[clone] = set(lung_vars)-set(pt_vars)

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


# Overview
fig, axs = plt.subplots(1,2,figsize=(6,4.5), gridspec_kw={'width_ratios': [7,3]})

df_ = df_clones.loc[:,~df_clones.columns.str.startswith('_')].reset_index(names='GBC').melt(id_vars='GBC')
order = ['PT', 'Lung', 'Enriched PT', 'Enriched lung', 'Enriched both', 'Acquired', 'Lost']
box(df_, x='variable', y='value', c='white', ax=axs[0], order=order)
strip(df_, x='variable', y='value', c='k', s=5, ax=axs[0], order=order)
format_ax(ax=axs[0], title='', rotx=45, ylabel='n MT-SNVs', reduced_spines=True)

df_ = df_clones.loc[:,df_clones.columns.str.startswith('_')].reset_index(names='GBC').melt(id_vars='GBC')
df_['variable'] = df_['variable'].map(lambda x: x[1:])
order = ['Enriched both', 'Acquired', 'Lost']
box(df_, x='variable', y='value', c='white', ax=axs[1], order=order)
strip(df_, x='variable', y='value', c='k', s=5, ax=axs[1], order=order)
format_ax(ax=axs[1], title='', rotx=45, ylabel='Mean AF +cells', reduced_spines=True)

fig.suptitle('~1 month of MT-SNVs clonal dynamics, in vivo')
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'MT_dynamics.png'), dpi=500)


##


# Enriched muts
muts = list(chain.from_iterable([ x for x in enriched_muts.values() ]))

X = np.zeros((len(enriched_muts)*2,len(muts)))
prevalence = np.zeros((len(enriched_muts)*2,len(muts)))
n_cells = np.zeros((len(enriched_muts)*2,len(muts)))
rows = []

i = 0
for clone in enriched_muts:
    cells_pt = afm.obs.query('GBC==@clone and Origin=="PT"').index
    X[i,:] = np.mean(afm[cells_pt,muts].X.A, axis=0)
    prevalence[i,:] = np.sum(afm[cells_pt,muts].layers['bin'].A>0, axis=0) / cells_pt.size
    n_cells[i,:] = np.sum(afm[cells_pt,muts].layers['bin'].A>0, axis=0)
    rows.append(f'{clone}_PT  ')
    i += 1
    cells_lung = afm.obs.query('GBC==@clone and Origin=="lung"').index
    X[i,:] = np.mean(afm[cells_lung,muts].X.A, axis=0)
    n_cells[i,:] = np.sum(afm[cells_lung,muts].layers['bin'].A>0, axis=0)
    prevalence[i,:] = np.sum(afm[cells_lung,muts].layers['bin'].A>0, axis=0) / cells_lung.size
    rows.append(f'{clone}_lung')
    i += 1

X = pd.DataFrame(X, index=rows, columns=muts)
prevalence = pd.DataFrame(prevalence, index=rows, columns=muts)
n_cells = pd.DataFrame(n_cells, index=rows, columns=muts)
n_cells = n_cells.loc[:,np.any(prevalence>0, axis=0)]
X = X.loc[:,np.any(prevalence>0, axis=0)]
prevalence = prevalence.loc[:,np.any(prevalence>0, axis=0)]

df_ = X.melt(var_name='mut', value_name='AF', ignore_index=False).reset_index(names='clone').merge(
    # n_cells.melt(var_name='mut', value_name='n_cells', ignore_index=False).reset_index(names='clone'),
    prevalence.melt(var_name='mut', value_name='prevalence', ignore_index=False).reset_index(names='clone'),
    on=['clone', 'mut']
)

# Viz
min_size = 0
max_size = 100  
from matplotlib.colors import Normalize
size_norm = Normalize(vmin=0, vmax=1) # 10

fig, ax = plt.subplots(figsize=(7,4))
sns.scatterplot(data=df_, x='mut', y='clone', hue='AF', size='prevalence', # size='n_cells', 
                legend='brief', ax=ax, sizes=(min_size, max_size), size_norm=size_norm, palette='mako_r', edgecolor='black')
format_ax(ax=ax, rotx=90)
ax = plt.gca()
ax.legend(
    bbox_to_anchor=(1.05, 1.05),  # Adjust the position
    loc='upper left',          # Location relative to the bounding box
    borderaxespad=0,           # Padding between the axes and legend box
    frameon=False
)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'enriched_MTs.png'), dpi=500) # 'acquired_MTs.png'), dpi=500)


##
