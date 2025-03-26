"""
Final AML MT-SNVs trees
"""

import os
import scanpy as sc
from mito_utils.phylo import *
from mito_utils.phylo_plots import *
matplotlib.use('macOSX')


##


# Set paths
sample = 'sAML1'
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro' 
path_results = os.path.join(path_main, f'results/AML/test_plot/{sample}')


##


# Colors
cmaps = {'tumor_tme_class_new' : {'tme':'b', 'malignant':'r'}}


##


# Open fig
fig, axs = plt.subplots(3,5,figsize=(15,9))

i = 0
for d,_,files in os.walk(path_results):
    if len(files)>0:

        # Tree
        with open(os.path.join(d, 'annotated_tree.pickle'), 'rb') as f:
            tree = pickle.load(f)
        metrics_df = pd.read_csv(os.path.join(d, 'tree_metrics.csv'), index_col=0).set_index('metric')
        afm = sc.read(os.path.join(d, 'afm.h5ad'))

        ax = axs.ravel()[i]
        plot_tree(tree, features=['tumor_tme_class_new'], ax=ax, categorical_cmaps=cmaps, colorstrip_width=20)
        n_cells = metrics_df.loc['n_cells','value']
        n_characters = metrics_df.loc['n_characters','value']
        median_support = metrics_df.loc['median_support','value']
        corr_distances = metrics_df.loc['corr_distances','value']
        title = f'cells x char: {n_cells} x {n_characters}\nSupp: {median_support:.2f}, Corr: {corr_distances:.2f}'
        ax.set(title=title)

        i+=1


# Save
fig.tight_layout()
plt.show()
fig.savefig(os.path.join(path_results, f'{sample}_trees.png'), dpi=300)


##