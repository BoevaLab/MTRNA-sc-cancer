# Quality control preprocessing script

This folder contains the script used to perform the quality control on the data used in the study.  

An example bash file to run the script is provided as `example-bash.sh`. Importantly, you must provide the name of the malignant cells and the TME cells in the cell type column. 

We also provide a list of the command lines we used to process the datasets in this study. 

**NOTE: the Chung et al. breast data used for Fig 1 and Suppl Fig S13 had to preprocessed separately because of the small amount of cells per patients.** We thus processed the full dataset rather than setting thresholds per patient. To do so, we added a "dummy_patient" column in the downloaded data and used that column as the `sample` column.

The filtered/processed data can also be provided directly to interested researchers. To do so, please contact josephine.yates@inf.ethz.ch or valentina.boeva@inf.ethz.ch. 

## Argument descriptions

- `data`: path to your .h5ad data. The data should be an h5ad object containing the raw counts or TPMs in .X and containing a column with the cell type and a column with the sample ID in .obs. Cell types will be further refined with CNVs. 
- `datatype`: whether the data contains raw counts or TPMs. For technologies like 10X, SeqWell, etc., use "10x"; for technologies like Smartseq-2 that do not only capture the 3' end and thus need gene-length correction, use "tpm". The "tpm" option will round the TPM counts to the nearest integer, and also skips additional hard thresholding filtering steps only adapted to 10X/SeqWell/etc. technologies.
- `savedir`: path to where you will save your data. If not already created, will be created.
- `gencodefile`: this file is used in the inferCNV step, as we need to have the information about gene position. We provide this file in the `auxiliary_data` folder of the repo.
- `samplecol`: the name of the column in your AnnData .obs where the sample ID is stored.
- `celltypecol`: the name of the column in your AnnData .obs where the cell type is stored.
- `malignantcells`: name of the malignant cells in your cell type column. If several, just write them with a space in between (e.g., "Malignant 1" "Malignant 2")
- `tmecells`: name of the TME cells in your cell type column. If several, just write them with a space in between (e.g., "TME 1" "TME 2")
- `nmintme`: this is the minimum number of TME cells present in a sample for the sample to be considered; otherwise, the sample will not be used. Setting this too low might cause issues with the inferCNV step. We used 50 as a default.
- `nmintmerefcat`: this is the minimul of TME cells of a specific type that need to exist in a sample for it to be considered a valid inferCNV category. We used 10 as a default.
- `leidenres`: this is the resolution used in the Leiden algorithm used at the CNV step. Larger values lead to more clusters. We used 3 as a default. 
