# MI_TO_analysis_repro
Reproducibility of all MI_TO paper downstream analyses. To reproduce the main results:

1. Clone and cd to this repository.
2. Install all the desired packages with mamba:

```bash
mamba env create -f envs/environment.yml -n MI_TO
```

3. Retrieve the project data folder at {zenodo link} and put it in the same folder of the cloned repo (i.e., path_main)
4. Create a results folder in path_main, and run the desired analysis from code (in progressive order, 1,2,3... and so on)

This repo is still in active development.
