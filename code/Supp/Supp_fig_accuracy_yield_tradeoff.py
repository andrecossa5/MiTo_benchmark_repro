"""
Supp Fig ... . Accuracy-yield trade-off.
Individual impact of hyper-parameters on metrics of interest.
"""

import os
import mito as mt
import matplotlib
import matplotlib.pyplot as plt
import plotting_utils as plu
from itertools import product
matplotlib.use('macOSX')


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_data = os.path.join(path_main, 'data', 'tune')
path_figures = os.path.join(path_main, 'results', 'figures', 'Supp')
path_results = os.path.join(path_main, 'results', 'others', 'Supp')


##



# Load tune results
df, metrics, options = mt.ut.format_tuning(path_data)

# Lists
metrics = ['ARI', 'NMI', 'AUPRC']
samples = ['MDA_clones', 'MDA_PT', 'MDA_lung']

# Params
plu.set_rcParams()

# Fig
fig, axs = plt.subplots(3,3,figsize=(10.2,8), sharex=True)

for i,(metric,sample) in enumerate(product(metrics, samples)):

    df_agg = (
        df
        .groupby(['sample', 'af_confident_detection'])
        .apply(lambda x: x[[metric, 'n_cells']]
        .median()).reset_index()
    )
    df_ = df_agg.query('sample==@sample')
    ax = axs.ravel()[i]
    ax.plot(df_['af_confident_detection'], df_[metric], 'bo--', label=metric)
    ax.set_ylabel(metric, color='b')
    ax.tick_params(axis='y', labelcolor='b')
    ax.set(xlabel='Min AF confident detection', title=sample)

    ax_ = ax.twinx()
    ax_.plot(df_['af_confident_detection'], df_['n_cells'], 'gx--', label='n_cells')
    if i>5:
        ax.set_xlabel('af_confident_detection')
    else:
        ax.set_xlabel('')
    ax_.set_ylabel('n_cells', color='g')
    ax_.tick_params(axis='y', labelcolor='g')

fig.tight_layout()
fig.savefig(os.path.join(path_figures, f'Supp_fig_accuracy_yield_tradeoff.pdf'))


##
