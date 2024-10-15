"""
Contribute of individual options to individual metrics.
"""

import os
from mito_utils.utils import *
from mito_utils.preprocessing import *
from mito_utils.plotting_base import *
matplotlib.use('macOSX')


##


# Get metrics
path_main = '/Users/IEO5505/Desktop/MI_TO/MI_TO_analysis_repro'
path_results = os.path.join(path_main, 'results', 'MI_TO_bench', 'benchmark')

# Main metrics and options df
df = pd.read_csv(os.path.join(path_results, 'main_benchmarking_df.csv'), index_col=0)


##


# 