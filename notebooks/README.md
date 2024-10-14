# Notebooks

This folder contains the notebooks used to perform the study. 
Placeholders were inserted as `/add/path/here/` in the files to replace paths. These will need to be replaced to run the analysis. 
Some of these placeholders indicate 

Any placeholder that starts with `/add/path/here/auxiliary_data` indicates a file that is provided in the `auxiliary_data` [folder](https://github.com/BoevaLab/MTRNA-sc-cancer/tree/main/auxiliary_data).
Replace these by wherever you saved the `auxiliary_data` folder.

In the following explanations, the `XXXX` replace names that are repeated across the notebook (usually the cancer type).

## `1. analyze_unfiltered_data.ipynb`
This notebook contains the analysis used to create Suppl Fig S2 and S12. 
- `adata = sc.read_h5ad("/add/path/here/XXXX.h5ad")`: path to where the data BEFORE preprocessing is saved.
- `filtered_adata = sc.read_h5ad("/add/path/here/filtered_data/XXXX/filtered_adata.h5ad")`: path to where the data AFTER preprocessing is saved.
- `fig.savefig("/add/path/here/figures/pre_vs_post/XXXX.svg", dpi=200, bbox_inches="tight")`:  path to where to save the figures.


## `2. celllevel_analysis.ipynb`
This notebook contains the analysis used to create Suppl Fig S3-S11.
- `full_resdir = pl.Path("/add/path/here/markers_highmt")`: where to save the results
- `adata = sc.read_h5ad("/add/path/here/filtered_data/XXXX/filtered_adata.h5ad")`: where the data after preprocessing data is saved.
- `fig.savefig("/add/path/here/figures/case_control/XXXX.png", dpi=200, bbox_inches="tight")`: where to save the figures

## `3. create_metacells.ipynb`
This notebook contains the analysis that computes the metacells, used in the rest of paper, notably Fig 1, Fig 3, Suppl Fig S17.
- `adata = sc.read_h5ad("/add/path/here/filtered_data/XXXX/filtered_adata.h5ad")`: where the data after preprocessing data is saved.
- `metacells.write("/add/path/here/metacell_data/XXXX/metacells.h5ad")`: where to save the metacells AnnData.
- `adata.obs.to_csv("/add/path/here/metacell_data/XXXX/adata_obs.csv")`: where to save the .obs of the AnnData linking metacells to the original cells.
- `egfr_df.to_csv("/add/path/here/egfr_expression_metacells/SCLC.csv")`: where to save the EGFR family gene expressions, useful for Suppl Fig S17.
- `diffs.sort_values("Diff. median").to_csv("/add/path/here/metabolic_dysregulation/mitocarta_met_res/SCLC_Chan_10X.csv")`: where to save the results of MitoCarta pathway dysregulation.
- `metacells.obs[['cleaned_celltype', 'Malignant', 'HighMT', 'CYP genes', 'UGT genes', 'GST genes', 'ABC transporters', 'O-Flanagan_dissociation_stress', 'O-Flanagan_dissociation_red', 'Van-den-Brink_dissociation_stress', 'Machado_dissociation_stress', 'Dissociation stress', 'Mito transfer']].to_csv("/add/path/here/info-metacells/sclc.csv")`: where the save the computed scores - these are useful to create plots.

## `4. bulk_vs_bulkified.ipynb`
This notebook contains the analysis for the bulk vs bulkified comparison, used in Fig 1 and Suppl Fig S13.

Data from Wu et al. can be found [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078).
Data from Chung et al. can be found [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi). Data was processed to separate the bulk measurements from the sc measurements.

- `bulk = pd.read_csv("/add/path/here/auxiliary_data/BRCA_Wu_fpkm.csv",index_col=0).T`: The FPKM transformed data from Wu et al.
- `adata = sc.read_h5ad("/add/path/here/filtered_data/Breast_Wu_10X/filtered_adata.h5ad")`: the path to the Wu et al. data AFTER preprocessing.
- `resdir = pl.Path(f"/add/path/here/results_bulk_vs_bulkified/XXXX/deg{deg}")`: where to save the results of the analysis.
- `mean_dir = pl.Path("/add/path/here/results_bulk_vs_bulkified/XXXX/")`: where the results were saved.
- `fig.savefig(f"/add/path/here/figures/bulk_vs_bulkified/XXXX_deg{deg}.svg", dpi=200, bbox_inches="tight")`: where to save the figures.
- `fig.savefig("/add/path/here/figures/bulk_vs_bulkified/XXXX_R2.svg", dpi=200, bbox_inches="tight")`: where to save the figures.
- `bulk = pd.read_csv("/add/path/here/auxiliary_data/Breast_Chung_bulk.csv",index_col=0)`: the FPKM transformed data from Chung et al.
- `adata = sc.read_h5ad("/add/path/here/filtered_data/Breast_Chung")`: the path to the Chung et al. data AFTER preprocessing. **NOTE:** The Chung et al. data is processed differently due to the small amount of cells per patients; see more info in the [preprocessing folder](https://github.com/BoevaLab/MTRNA-sc-cancer/tree/main/preprocessing). 


