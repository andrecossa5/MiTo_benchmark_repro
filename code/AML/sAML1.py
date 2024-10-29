"""
sAML1 analysis
"""

import os
from mito_utils.utils import *
from mito_utils.clustering import *
from mito_utils.preprocessing import *
from mito_utils.dimred import *
from mito_utils.plotting_base import *
from mito_utils.diagnostic_plots import *
from mito_utils.embeddings_plots import *
from mito_utils.phylo_plots import *
matplotlib.use('macOSX')


##


# Args
sample ='sAML1'
job_id = 'tuningf317e5900d'
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro' 
path_tuning = os.path.join(path_main, f'results/AML/tuning/{sample}')
path_afm = os.path.join(path_main, f'data/AML/AFMs/{sample}/afm.h5ad')
path_results = os.path.join(path_main, f'results/AML/phylo_inference/{sample}')
# path_coverage = os.path.join(path_main, f'data/MI_TO_bench/AFMs/mito_preprocessing/{sample}/coverage.txt.gz')
# wobble_table = os.path.join(path_main, 'data/MI_TO_bench/miscellanea/formatted_table_wobble.csv')
# weng_2024_ref = os.path.join(path_main, 'data/MI_TO_bench/miscellanea/weng2024_mut_spectrum_ref.csv')


afm = sc.read(path_afm)
with open(os.path.join(path_tuning, f'{job_id}_stats.pickle'), 'rb') as f:
    d = pickle.load(f)

afm = filter_cells(afm, cell_filter='filter2')

afm.obs['lymph_vs_myelo'] = np.select(
    [afm.obs['aggregated_ct'].isin(['T_CD8', 'T_CD4']), 
     afm.obs['aggregated_ct'].isin(['Early_myeloid', 'HSC_MPP', 'B_early'])],
    ['lympho', 'progenitors'], default='Mature myelo'
) 
afm = filter_afm(afm, filtering=None, 
    cells=d['cells'].to_list(), 
    variants=d['vars'].to_list(), 
    bin_method='vanilla', 
    binarization_kwargs={'t_vanilla':0, 'min_AD':d['genotyping']['min_AD']}
)
tree = build_tree(afm, metric='jaccard', bin_method='vanilla', solver='NJ')
_, mut_nodes, mut_order = cut_and_annotate_tree(tree)


##


tree.collapse_mutationless_edges(True)


##

# fig, axs = plt.subplots(1,2,figsize=(10,5))
# 
# cmaps = {
#     'tumor_tme_class_new' : {'malignant':'red', 'tme':'b', 'undefined':'grey'},
#     'lymph_vs_myelo' : { 'lympho' : 'green', 'progenitors' : 'orange', 'Mature myelo' : 'darkred'},
#     'MT_clone' : create_palette(tree.cell_meta, 'MT_clone', ten_godisnot)
# }
# 
# ax = axs[0]
# plot_tree(tree, features=['tumor_tme_class_new', 'lymph_vs_myelo', 'MT_clone'], 
#           categorical_cmaps=cmaps, ax=ax, colorstrip_width=1)
# ax.set(title=f'n_cells: {afm.shape[0]}, n_vars {afm.shape[1]}')
# 
# ax = axs[1]


fig, ax = plt.subplots(figsize=(5,5))
plot_tree(tree, features=mut_order, layer='transformed', ax=ax, colorstrip_width=2, 
          colorstrip_spacing=0, orient='down')
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'{job_id}_trees.png'), dpi=500)

# plt.show()








# fig.savefig(os.path.join(path_results, sample, f'{job_id}_tree.png'))


fig, ax = plt.subplots(figsize=(4.8,5))
order=leaves_list(linkage(afm.obsp['distances'].A))
ax.imshow(afm.obsp['distances'].A[np.ix_(order, order)], cmap='inferno_r')
fig.tight_layout()
plt.show()

# fig.savefig(os.path.join(path_results, sample, f'{job_id}_distances.png'), dpi=500)


fig, ax = plt.subplots(figsize=(4,4))
ax.plot(
    np.sum(afm[afm.obs['tumor_tme_class_new']=='tme'].layers['bin'].A==1, axis=0) / np.sum(afm.obs['tumor_tme_class_new']=='tme'),
    np.sum(afm[afm.obs['tumor_tme_class_new']=='malignant'].layers['bin'].A==1, axis=0) / np.sum(afm.obs['tumor_tme_class_new']=='malignant'),
    'ko'
)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set(xlim=(-0.1,1.1), ylim=(-0.1,1.1))
format_ax(ax=ax, xlabel='TME', ylabel='Tumor')
fig.tight_layout()
fig.savefig(os.path.join(path_results, sample, f'{job_id}_tme_tumor.png'))
plt.close()

