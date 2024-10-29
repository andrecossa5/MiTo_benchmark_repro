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
path_data = os.path.join(path_main, 'results', 'MI_TO_bench', 'benchmark', 'NJ')
path_results = os.path.join(path_main, 'results', 'MI_TO_bench', 'benchmark')


# Read
sample = 'MDA_lung'

n_positive = []
mean_AD_in_positives = []
variants = []
for x in os.listdir(os.path.join(path_data, sample)):
    afm = sc.read(os.path.join(path_data, sample, x, 'afm.h5ad'))
    variants += afm.var_names.to_list()
    n_positive += afm.var['Variant_CellN'].to_list()
    mean_AD_in_positives += afm.var['mean_AD_in_positives'].to_list()

df = pd.DataFrame({'n_positive':n_positive, 'mean_AD_in_positives':mean_AD_in_positives}, index=variants).drop_duplicates().reset_index()


##


fig, ax = plt.subplots(figsize=(4.5,4.5))

ax.set_yscale('log', base=2)
ax.set_xscale('log', base=2)

sns.kdeplot(data=df, x='n_positive', y='mean_AD_in_positives', fill=False, ax=ax)
ax.plot(df['n_positive'], df['mean_AD_in_positives'], marker='o', linestyle='', color='darkorange', markersize=5, markeredgecolor='k')
xticks = [0,1,2,5,10,20,40,80,200,500]
yticks = [0,1,2,4,8,16,32,64,264]
ax.xaxis.set_major_locator(FixedLocator(xticks))
ax.yaxis.set_major_locator(FixedLocator(yticks))

def integer_formatter(val, pos):
    return f'{int(val)}'

ax.xaxis.set_major_formatter(FuncFormatter(integer_formatter))
ax.yaxis.set_major_formatter(FuncFormatter(integer_formatter))

ax.set(title=f'{sample}: {df.shape[0]} MT-SNVs', xlabel='n +cells', ylabel='Min n ALT in +cells')

fig.tight_layout()
fig.savefig(os.path.join(path_results, f'{sample}_sweet_spot.png'), dpi=500)


##