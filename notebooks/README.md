# Notebooks

This folder contains the notebooks used to perform the study. 
Placeholders were inserted as `/add/path/here/` in the files to replace paths. These will need to be replaced to run the analysis. 
Some of these placeholders indicate 

Any placeholder that starts with `/add/path/here/auxiliary_data` indicates a file that is provided in the `auxiliary_data` [folder](https://github.com/BoevaLab/MTRNA-sc-cancer/tree/main/auxiliary_data).

In the following explanations, the `XXXX` replace names that are repeated across the notebook (usually the cancer type).

## `analyze_unfiltered_data.ipynb`
This notebook contains the analysis used to create Suppl Fig S2 and S12. 
- `adata = sc.read_h5ad("/add/path/here/XXXX.h5ad")`: path to where the data BEFORE preprocessing is saved.
- `filtered_adata = sc.read_h5ad("/add/path/here/filtered_data/XXXX/filtered_adata.h5ad")`: path to where the data AFTER preprocessing is saved.
- `fig.savefig("/add/path/here/figures/pre_vs_post/XXXX.svg", dpi=200, bbox_inches="tight")`:  path to where to save the figures.


## `celllevel_analysis.ipynb`
This notebook contains the analysis used to create Suppl Fig S3-S11.
- `full_resdir = pl.Path("/add/path/here/markers_highmt")`: where to save the results
- `adata = sc.read_h5ad("/add/path/here/filtered_data/XXXX/filtered_adata.h5ad")`: where the data after preprocessing data is saved.
- `fig.savefig("/add/path/here/figures/case_control/XXXXX.png", dpi=200, bbox_inches="tight")`: where to save the figures

##




