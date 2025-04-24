"""
Bench cell_filter
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


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
# path_figures = os.path.join(path_main, 'results', 'figures', 'Fig2')
# path_results = os.path.join(path_main, 'results', 'others', 'Fig2')



# 1. Per site coverage metrics, filter1 and filter2 -------------------------- 

# Path data
path_data = os.path.join(path_main, 'data', 'general', 'AFMs', 'maegatk')

# Params
plu.set_rcParams()
matplotlib.rcParams.update({'figure.dpi':150})

# Read afm
sample = 'MDA_PT'
afm = sc.read(os.path.join(path_data, f'afm_{sample}.h5ad'))

# Add tests columns
median_target_site_coverage = 25
frac_target_site_covered = .75
mean_cov_all = 20
nmads = 5

# No filter
afm.obs['No filter'] = 'Filtered'

# filter1
x = afm.obs['mean_site_coverage']       
median = np.median(x)
MAD = np.median(np.abs(x-median))
test = (x>=mean_cov_all) & (x<=median+nmads*MAD)
afm.obs['filter1'] = 'Discarded'
afm.obs.loc[test, 'filter1'] = 'Filtered'

# filter2
test1 = afm.obs['median_target_site_coverage'] >= median_target_site_coverage
test2 = afm.obs['frac_target_site_covered'] >= frac_target_site_covered
afm.obs['filter2'] = 'Discarded'
afm.obs.loc[(test1) & (test2), 'filter2'] = 'Filtered'


##


fig, axs = plt.subplots(1,3,figsize=(7.5,3),sharey=True)

cell_filters = ['No filter', 'filter1', 'filter2']
cmap = {'Discarded':'grey', 'Filtered':'darkred'}

for i,cell_filter in enumerate(cell_filters):

    ax = axs.ravel()[i]
    plu.scatter(afm.obs, 
                x='median_target_site_coverage', 
                y='frac_target_site_covered',
                by=cell_filter, categorical_cmap=cmap,
                size=10, alpha=.3, ax=ax)
    
    x_median = afm.obs.loc[afm.obs[cell_filter]=='Filtered', 'median_target_site_coverage'].median()
    y_median = afm.obs.loc[afm.obs[cell_filter]=='Filtered', 'frac_target_site_covered'].median()
    n_cells = np.sum(afm.obs[cell_filter]=='Filtered')

    plu.format_ax(ax=ax, 
                  title=f'Cell filter: {cell_filter}',
                  xlabel='MT-target site coverage',
                  ylabel='% MT-target site covered' if i==0 else '', 
                  reduced_spines=True)

    ax.axvline(x=x_median, linewidth=.5, c='darkred', linestyle='--')
    ax.axhline(y=y_median, linewidth=.5, c='darkred', linestyle='--')
    ax.text(.25, .55, f'MT-target site coverage: {x_median:.2f}', transform=ax.transAxes, fontsize=7)
    ax.text(.25, .50, f'% MT-target site covered: {y_median:.2f}', transform=ax.transAxes, fontsize=7)
    ax.text(.25, .45, f'n cells: {n_cells:.2f}', transform=ax.transAxes, fontsize=7)

plu.add_legend(cmap, ax=axs[0], loc='lower right', label='Cell', bbox_to_anchor=(1,0))

fig.tight_layout()
plt.show()


##


# Last. 
L = []
afm.obs['GBC'] = afm.obs['GBC'].astype(str)

for cell_filter in cell_filters:
    for min_cell_number in [0,5]:
        n_clones = (
            afm.obs.loc[afm.obs[cell_filter]=='Filtered', 'GBC']
            .value_counts()
            .loc[lambda x: x>=min_cell_number]
            .index.size
        )
        L.append([cell_filter, min_cell_number, n_clones])

df_plot = pd.DataFrame(L, columns=['cell_filter', 'min_cell_number', 'n clones'])
df_plot['min_cell_number'] = df_plot['min_cell_number'].astype(str)
by_order = ['No filter', 'filter1', 'filter2']
cmap = plu.create_palette(df_plot, 'cell_filter', order=by_order[1:], palette='Blues')
cmap = {**{'No filter':'lightgrey'},**cmap}

fig, ax = plt.subplots(figsize=(4,3))
plu.bar(df_plot, x='min_cell_number', y='n clones', 
        by='cell_filter', ax=ax, categorical_cmap=cmap, by_order=by_order)
plu.format_ax(ax=ax, ylabel='n clones', xlabel='Min cell number', reduced_spines=True)
plu.add_legend(cmap, label='Cell filter', ax=ax)

fig.subplots_adjust(right=.6, bottom=.15, left=.15, top=.85)
plt.show()




##


# 2. Extended summary  -------------------------- 

path_data = os.path.join(path_main, 'data', 'bench', 'tune_cell_filter')


L = []
for folder,_,files in os.walk(path_data):
    if any([ x.startswith('all') for x in files]):
        df,metrics,options = mt.ut.format_tuning(folder)
        L.append(df)
df = pd.concat(L)
df['cell_filter'].loc[lambda x: x.isna()] = 'No filter'

# Define options and metrics
varying_options = (df[options+['cell_filter']].nunique()).loc[lambda x:x>1].index.to_list()
metrics_of_interest = ['ARI', 'NMI', 'corr', 'AUPRC', 'n_cells', 'n_vars', 'n_GBC_groups']

(
    df.groupby(['sample']+varying_options)
    [['ARI', 'NMI', 'corr', 'AUPRC', 'n_cells', 'n_GBC_groups']]
    .median()
)


##


# Cell filters 
fig, axs =  plt.subplots(2,len(metrics_of_interest),figsize=(14,5.5))

x_order = ['MDA_clones', 'MDA_lung', 'MDA_PT']
by_order = ['No filter', 'filter1', 'filter2']

cmap = plu.create_palette(df, 'cell_filter', order=by_order[1:], palette='Blues')
cmap = {**{'No filter':'lightgrey'},**cmap}

for i,metric in enumerate(metrics_of_interest):
    plu.bar(
        df.query('min_cell_number=="0"'), 
        x='sample', y=metric, 
        by='cell_filter',
        x_order=x_order,
        by_order=by_order,
        categorical_cmap=cmap,
        ax=axs[0,i]
    )
    plu.format_ax(ax=axs[0,i], xticks=['' for _ in range(3) ], 
                  xlabel='', ylabel=metric, reduced_spines=True)

for i,metric in enumerate(metrics_of_interest):
    plu.bar(
        df.query('min_cell_number=="5"'), 
        x='sample', y=metric, 
        by='cell_filter',
        x_order=x_order,
        by_order=by_order,
        categorical_cmap=cmap,
        ax=axs[1,i]
    )
    plu.format_ax(ax=axs[1,i], rotx=90, xlabel='', ylabel=metric, reduced_spines=True)

fig.tight_layout()
plu.add_legend(cmap, label='Cell filter', ax=axs[0,6])
fig.subplots_adjust(top=.90, bottom=.2, right=.78, left=.1, wspace=.9)
plt.show()


##




# fig, axs = plt.subplots(1,len(metrics_of_interest),figsize=(14,4))

# x_order = ['0', '5']
# by_order = ['No filter', 'filter1', 'filter2']

# cmap = plu.create_palette(df, 'cell_filter', order=by_order[1:], palette='Blues')
# cmap = {**{'No filter':'lightgrey'},**cmap}

# for i,metric in enumerate(metrics_of_interest):
#     plu.strip(df, x='min_cell_number', y=metric, 
#             x_order=x_order, by_order=by_order,
#             by='cell_filter', categorical_cmap=cmap, ax=axs[i])
#     kwargs = {   
#         'boxprops' : {'linewidth': 0}, 
#         'medianprops': {"color": "black", "linewidth": 1.5},
#         'whiskerprops':{'linewidth':0},
#         'fill':False
#     }
#     plu.box(df, x='min_cell_number', y=metric,
#             x_order=x_order, by_order=by_order, kwargs=kwargs,
#             by='cell_filter', ax=axs[i])#, categorical_cmap=cmap)
#     plu.format_ax(ax=axs[i], xlabel='n cells', title=metric, reduced_spines=True)
#     if i==0:
#         plu.format_ax(ax=axs[i], ylabel='Value')#, xlabel='n cells', title=metric, reduced_spines=True)
# 
# plu.add_legend(cmap, label='Cell filter', ax=axs[i])
# fig.subplots_adjust(top=.85, bottom=.15, right=.78, left=.1, wspace=.8)
# plt.show()


##