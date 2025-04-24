"""
Supp Fig ... . Meta-analysis MiTo benchmark metrics and options
Individual impact of hyper-parameters on metrics of interest.
"""

import os
import pandas as pd
import mito as mt
import matplotlib
import matplotlib.pyplot as plt
import plotting_utils as plu
from lightgbm import LGBMRegressor
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



##


# Prep data
df_options = df.loc[:,options].copy()
variable_options = df_options.nunique().loc[lambda x: x>1].index

X = df_options.loc[:,variable_options]
for x in X.columns:
    try:
        X[x] = X[x].astype(float)
    except:
        X[x] = X[x].astype('category')

##

# Select test metrics
metric_tests = [
    'AUPRC', 'ARI', 'NMI', 'freq_lineage_biased_muts', 'mean_CI', 'corr',
    'n_cells', 'median_n_vars_per_cell'
]
names = {
    'freq_lineage_biased_muts' : "% lineage-biased \n MT-SNVs",
    'n_cells' : 'n cells',
    'median_n_vars_per_cell' : 'n MT-SNVs per cell',
    'mean_CI' : 'CI',
    'corr' : 'Tree-char dists correlation'
}

# Here we go
plu.set_rcParams()

fig, axs = plt.subplots(2,4,figsize=(14,6))

for i, metric in enumerate(metric_tests):

    ax = axs.ravel()[i]
    y = df.loc[df_options.index, metric]

    model = LGBMRegressor(n_estimators=100, learning_rate=0.01, 
                          min_data_in_leaf=1, importance_type='split')
    model.fit(X, y)
    
    feature_importances = model.feature_importances_
    df_ = (
        pd.Series(feature_importances, X.columns)
        .sort_values(ascending=False)
        .to_frame('Importance')
    )

    ax.hlines(y=df_.index, xmin=0, xmax=df_['Importance'], color='k', linewidth=1.5, alpha=.7)
    ax.plot(df_['Importance'], df_.index, "o", color='darkred', markersize=8, markeredgecolor='k')
    ax.invert_yaxis()
    title = metric if not metric in names else names[metric]
    plu.format_ax(ax=ax, title=title, reduced_spines=True, xlabel='Feature importance')

fig.tight_layout()
fig.savefig(os.path.join(path_figures, 'Supp_fig_feature_importance_meta_analysis.pdf'))


## 