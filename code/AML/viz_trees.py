"""
Visualizing final AML MT-SNVs spaces and trees.
"""

import os
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.dimred import *
from mito_utils.plotting_base import *
from mito_utils.diagnostic_plots import *
from mito_utils.embeddings_plots import *
from mito_utils.phylo import *
from mito_utils.phylo_plots import *
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro' 
path_results = os.path.join(path_main, f'results/AML/final_trees/')
weng_2024_ref = os.path.join(path_main, 'data/MI_TO_bench/miscellanea/weng2024_mut_spectrum_ref.csv')


##


# Collect all plotting combos
# sample ='sAML1'
# job = 'tuning7b619458a9'
L = []
for d,_,file in os.walk(path_results):
    if len(file)>0 and file[0].startswith('ann'):
        sample,job = d.split('/')[-2:]
        L.append([sample,job])

# Here we go
for sample,job in L:

    # Load tree
    path_results = os.path.join(path_main, f'results/AML/final_trees/{sample}/{job}')
    with open(os.path.join(path_results, 'annotated_tree.pickle'), 'rb') as f:
        tree = pickle.load(f)

    # Mut profile
    ref_df = pd.read_csv(weng_2024_ref, index_col=0)
    fig = mut_profile(tree.layers['raw'].columns, ref_df=ref_df, figsize=(5,2.5))
    fig.tight_layout()
    fig.savefig(os.path.join(path_results, 'mut_profile.png'), dpi=500)


    ##


    # Annotate to get mut_nodes and mutation_order
    _, mut_nodes, mutation_order = MiToTreeAnnotator(tree)
    cmaps = {
        'tumor_tme_class_new' : {'tme':'b', 'malignant':'r'},
        'MT_clone' : create_palette(tree.cell_meta, 'MT_clone', sc.pl.palettes.default_102)
    }


    ##


    # Optional
    # tree.collapse_mutationless_edges(True)
    # calculate_corr_distances(tree)


    ##

    fig, axs = plt.subplots(1,2,figsize=(8,4.5))

    ax = axs[0]
    plot_tree(tree, ax=ax, feature_internal_nodes='support', internal_node_subset=mut_nodes,
              cmap_internal_nodes=('Spectral_r',50,80), internal_node_kwargs={'markersize':5})
    # corr_distances = tree_metrics.loc["corr_distances", 'value']
    # supp_mut_nodes = np.median(get_supports(tree, subset=mut_nodes))
    # ax.set(title=f'Mut-assigned nodes support: {supp_mut_nodes}\nTree-MT-SNVs distances corr: {corr_distances:.2f}')
    # add_cbar(get_supports(tree), ax=ax, palette='Spectral_r', vmin=50, vmax=80, label='Support', layout='h3')

    ax = axs[1]
    plot_tree(tree, features=['tumor_tme_class_new', 'MT_clone'], ax=ax, categorical_cmaps=cmaps, colorstrip_width=3)
    # ARI = tree_metrics.loc["ARI", 'value']
    # NMI = tree_metrics.loc["NMI", 'value']
    # ax.set(title=f'ARI: {ARI:.2f}, NMI: {NMI:.2f}')

    fig.tight_layout()
    fig.savefig(os.path.join(path_results, 'phylo_main.png'), dpi=500)


##


# Main mut trees
# fig, axs = plt.subplots(1,2,figsize=(15,8), gridspec_kw={'wspace': 0.4})
# 
# plot_tree(tree, ax=axs[0], 
#     colorstrip_spacing=.000001, colorstrip_width=2,
#     orient='down',
#     features=['GBC', 'MT_clone']+mutation_order, layer='raw',
#     categorical_cmaps=cmaps,
#     feature_label_size=10, feature_label_offset=2,
# )
# add_cbar(tree.layers['transformed'].values.flatten(), palette='mako', label='AF', ticks_size=8, label_size=9,
#          vmin=.01, vmax=.1,
#          ax=axs[0], layout=( (1.02,.3,.02,.2), 'right', 'vertical' ))
# 
# plot_tree(tree, ax=axs[1],
#     colorstrip_spacing=.000001, colorstrip_width=2,
#     orient='down',
#     features=['GBC', 'MT_clone']+mutation_order, layer='transformed',
#     feature_label_size=10, feature_label_offset=2,
#     categorical_cmaps=cmaps,
#     internal_node_subset=mut_nodes,
#     internal_node_kwargs={'markersize':5, 'c':'darkred'}, show_internal=True
# )
# add_legend(label='Genotype', ax=axs[1], colors={'REF':'b', 'ALT':'r'}, loc='center left', bbox_to_anchor=(1,.4),
#            ticks_size=8, artists_size=10, label_size=9)
# 
# fig.tight_layout()
# fig.savefig(os.path.join(path_results, 'mut_tree.png'))


##