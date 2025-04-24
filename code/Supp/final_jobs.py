"""
Select final jobs for clonal inference benchmark.
"""

import os
import pandas as pd
import mito as mt


##


# Set paths
path_main = '/Users/IEO5505/Desktop/MI_TO/MiTo_benchmark_repro'
path_results = os.path.join(path_main, 'results', 'others', 'Fig2')


##


# Format tune output
path_data = os.path.join(path_main, 'data', 'tune')
df, metrics, options = mt.ut.format_tuning(path_data)
df['frac_unassigned'] = df['unassigned'] / df['n_cells']

##

# Filters
d = {
    'MDA_clones' : {'n_cells':250, 'n_GBC_groups':6},
    'MDA_PT' : {'n_cells':1000, 'n_GBC_groups':30},
    'MDA_lung' : {'n_cells':1000, 'n_GBC_groups':10},
}
frac_unassigned = .1

# Here we go
path_ = '/data/cossa2_dare/MI_TO_bench/data/AFMs/maegatk'
samples = ['MDA_clones', 'MDA_PT', 'MDA_lung']

for sample in samples:
    df_selected, df_final = mt.ut.select_jobs(
        df, sample, d[sample]['n_cells'], d[sample]['n_GBC_groups'], frac_unassigned
    )
    df_selected.to_csv(os.path.join(path_results, f'{sample}_filtered_jobs.csv'))
    (   
        df_final[['job_id']]    
        .assign(sample=sample, afm=os.path.join(path_, sample, 'afm.h5ad')) 
        .to_csv(os.path.join(path_results, f'{sample}_final_jobs.csv'), index=False) 
    )


##