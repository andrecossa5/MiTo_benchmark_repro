"""
Find sweet spot in MT-SNVs.
"""

import os
import pandas as pd
from mito_utils.preprocessing import *
from mito_utils.plotting_base import *
from mito_utils.diagnostic_plots import *
matplotlib.use('macOSX')


##


# Set paths and code
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro'
path_data = os.path.join(path_main, 'results', 'MI_TO_bench', 'phylo_inference', 'trees')
path_results = os.path.join(path_main, 'results', 'MI_TO_bench', 'phylo_inference')

os.listdir(os.path.join(path_data))

# Here we go
fig, axs = plt.subplots(1,3,figsize=(11,4))

for i,sample in enumerate(['MDA_clones', 'MDA_PT', 'MDA_lung']):
    
    n_positive = []
    mean_AD_in_positives = []
    variants = []

    for job in os.listdir(os.path.join(path_data, sample)):
        afm = sc.read(os.path.join(path_data, sample, job, 'afm.h5ad'))
        variants += afm.var_names.to_list()
        n_positive += (afm.layers['bin'].A==1).sum(axis=0).tolist()
        mean_AD_in_positives += np.nanmean(np.where(afm.layers['bin'].A==1, afm.layers['AD'].A, np.nan), axis=0).tolist()

    df = (
        pd.DataFrame(
            {'n_positive':n_positive, 'mean_AD_in_positives':mean_AD_in_positives}, 
            index=variants
        )
        .reset_index(names='var')
    )
    df.to_csv(os.path.join(path_results, f'{sample}_top_10_var_subset.csv'))

    ax = axs.ravel()[i]
    sns.kdeplot(data=df, x='n_positive', y='mean_AD_in_positives', fill=False, color='#41767F', ax=ax)
    sns.kdeplot(data=df, x='n_positive', y='mean_AD_in_positives', fill=True, color='#41767F', alpha=.7, ax=ax)
    n_vars = df['var'].unique().size
    median_n_positives = df['n_positive'].median()
    median_mean_AD_in_positives = df['mean_AD_in_positives'].median()
    format_ax(title=f'{sample}: {n_vars} MT-SNVs\nn +cells: {median_n_positives:.2f}, mean ALT +cells: {median_mean_AD_in_positives:.2f}', xlabel='n +cells', ylabel='Mean n ALT UMIs in +cells', ax=ax, reduced_spines=True)
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    ax.vlines(x=median_n_positives, ymin=y_min, ymax=median_mean_AD_in_positives, colors='red', linestyles='dashed')
    ax.hlines(y=median_mean_AD_in_positives, xmin=x_min, xmax=median_n_positives, colors='red', linestyles='dashed')


    ax.plot(median_n_positives, median_mean_AD_in_positives, 'rx', markersize=10)


fig.tight_layout()
fig.savefig(os.path.join(path_results, f'{sample}_sweet_spot.pdf'))


##