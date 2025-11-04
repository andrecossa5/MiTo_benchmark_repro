# MiTo_benchmark_repro
This repo holds the code to reproduce all [MiTo](https://www.biorxiv.org/content/10.1101/2025.06.17.660165v1) paper downstream analyses and figures (main and supplementary).

## Collect source data
To reproduce the paper analyses, start from downloading source data from [Zenodo](https://zenodo.org/records/17425550).

```bash
wget https://zenodo.org/records/17425550/files/source_data.tar.gz
tar -xzf source_data.tar.gz
tree source_data
```

this folder stores all preprocessed `data` and the associated `results` folders supporting the manuscript findings.

## Reproduce environment
To run this code, we need to install `MiTo` and its dependencies. See [MiTo](https://github.com/andrecossa5/MiTo).

## Reproduce analyses
The `code` directory stores all the code to reproduce the main (Fig. 2-5) and Supplementary (Supp 2-20) figures. To acess that, clone the repo:

```bash
git clone 
cd MiTo_benchmark_repro
tree code
```

The `code` folder structure should appear like this:

```
code
├── Fig2
│   └── Fig_2.py
├── Fig3
│   └── Fig_3.py
├── Fig4
│   ├── 1_gene_expression.py
│   ├── 2_Fig_4bc.py
│   ├── 3_Fig_4de.py
│   ├── 4_Fig_4f-l.py
│   └── 5_Fig_4mn.py
├── Fig5
│   └── Fig_5.py
└── Supp
    ├── Supp_Fig_10.py
    ├── Supp_Fig_11.py
    ├── Supp_Fig_12.py
    ├── Supp_Fig_13.py
    ├── Supp_Fig_14.py
    ├── Supp_Fig_15.py
    ├── Supp_Fig_16.py
    ├── Supp_Fig_17.py
    ├── Supp_Fig_18.py
    ├── Supp_Fig_19.py
    ├── Supp_Fig_2.py
    ├── Supp_Fig_20
    │   ├── 1_preprocess_Cas9_GEX_HASH.py
    │   ├── 2_make_AFM_Cas9.py
    │   ├── 3_redeem_data_cleaning.R
    │   ├── 4_make_AFM_MiTo.py
    │   ├── 5_Supp_Fig_20a.py
    │   ├── 6_Supp_Fig_20b.py
    │   ├── 7_Supp_Fig_20c.py
    │   └── 8_Supp_Fig_20d.py
    ├── Supp_Fig_3.py
    ├── Supp_Fig_5.py
    ├── Supp_Fig_6.py
    ├── Supp_Fig_7.py
    ├── Supp_Fig_8.py
    └── Supp_Fig_9.py
```

Each figure has its own script (or subdirectory), containing one (or more) .py script that read (write) from (to) the downloaded `source data` directory. To run the scripts on your machine, please set the `path_main` variable to match your `source data` directory _absolute_ path.

Each figure can be reproduced independently from the others. In case of more complex figures (i.e., Fig.4 and Supp. Fig. 20), to reproduce individual panels we must go through the corresponding scripts in order (e.g., for Fig. 4: `1_gene_expression.py`, `2_Fig_4bc.py`, `3_Fig_4de.py`, etc.)








