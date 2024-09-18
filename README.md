# MI_TO_analysis_repro
Reproducibility of all MI_TO paper downstream analyses. To reproduce the main results:

1. Clone and cd to this repository.
2. Install all the desired packages with conda/mamba:

```bash
mamba env create -f envs/environment.yml -n MI_TO
```

3. Manually install 3 additional packages:

```bash
pip install git+https://github.com/YosefLab/Cassiopeia@master#egg=cassiopeia-lineage
pip install -U git+https://github.com/single-cell-genetics/MQuad
# pip install mito_utils==1.0.0                                                         # Still under deveopment
# For now, clone git@github.com:andrecossa5/mito_utils.git and cd to it
# then: mamba develop .
```

4. Verify correct environment configuration:


```bash
mamba activate MI_TO
```
```python
import mito_utils
```

4. Retrieve the project data folder at {zenodo link} and put it in the same folder of the cloned repo (i.e., path_main)
5. Create a results folder in path_main, and run the desired analysis from code (in progressive order, 1,2,3... and so on)

N.B.: this repo is still in active development.
