# MTRNA-sc-cancer
Repository containing all the code necessary for reproducing the paper "Filtering cells  with high mitochondrial content removes viable metabolically altered malignant cell populations in cancer single-cell studies" (Yates, Kraft, and Boeva)

## How to reproduce the analysis

First, you will have to create a conda environment with the correct requirements. 
You will find a YAML file with the environment used to run the analysis. 

If you do not have Anaconda, you can download it here. 
You can create then a conda environment from the file using 
```
conda env create --name mtrna-env --file mtrna-env.yml
```
Activate the conda environment with 
```
conda activate mtrna-env
```


The first step of the analysis consists in downloading the data from the original source and transforming it so that it is saved as an .h5ad file to preprocess it. Details are given in [preprocessing](https://github.com/BoevaLab/MTRNA-sc-cancer/tree/main/preprocessing). 

Then, you can run the notebooks in the order indicated. Placeholders must be replaced in the files - description of the placeholders can be found in [notebooks](https://github.com/BoevaLab/MTRNA-sc-cancer/tree/main/notebooks).

## Where to download raw data

The data used in the paper can be found in the following places:
The single-cell studies used in this study can be downloaded from: 
- The Gene Expression Omnibus (GEO) website: Breast cancer, Wu et al., at GSE176078 and EGAD00001007495; Pancreatic ductal adenocarcinoma, Steele et al., at GSE155698; Prostate cancer, Song et al., at GSE176031; Nasopharyngeal carcinoma, Chen et al., at GSE150430; Breast cancer, Chung et al., at GSE75688
- The Broad single-cell portal: Metastatic Pancreatic cancer, Raghavan et al. at https://singlecell.broadinstitute.org/single_cell/study/SCP1644/microenvironment-drives-cell-state-plasticity-and-drug-response-in-pancreatic-cancer; Renal clear cell cancer, Bi et al. at https://singlecell.broadinstitute.org/single_cell/study/SCP1288/tumor-and-immune-reprogramming-during-immunotherapy-in-advanced-renal-cell-carcinoma#study-summary
- The cancer cell atlas (3CA): Small cell lung cancer, Chan et al., at https://www.weizmann.ac.il/sites/3CA/lung; Lung adenocarcinoma, Bischoff et al., https://www.weizmann.ac.il/sites/3CA/lung; Uveal Melanoma, Durante et al., https://www.weizmann.ac.il/sites/3CA/othermodels
Zenodo: mtDNA-linked single-cell, Kim et al., https://doi.org/10.5281/zenodo.10498240
- Single-cell lineage tracking data: Kuramochi cell line at GSE223003; MDAMB468 cell line at GSE228382

The bulk data used in this study can be downloaded from:
- The Gene Expression Omnibus (GEO) website: Breast cancer, Wu et al., at GSE176078; Breast cancer, Chung et al., at GSE75688
- The Cancer Cell Line Encyclopedia (CCLE): for the CCLE RNA-seq and drug sensitivity data https://depmap.org/portal/data_page/?tab=allData

The two samples processed with spatial transcriptomics method Visium HD are freely available on the 10X website: DCIS https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-human-breast-cancer-fresh-frozen and LUAD at https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-human-lung-cancer-post-xenium-expt.  
