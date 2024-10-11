import pandas as pd
import numpy as np
import scanpy as sc

import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import median_abs_deviation

from tqdm.notebook import tqdm

import infercnvpy as cnvpy

import os
import pathlib
import argparse
import logging
from typing import Optional, List

LOGGER = logging.getLogger(__name__)


def create_parser() -> argparse.ArgumentParser:
    """Creates the CLI parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "data", type=pathlib.Path, help="HDF5 file containing the dataset."
    )
    parser.add_argument(
        "--datatype",
        type=str,
        choices=["tpm", "10x"],
        help="path to where the data will be saved",
    )
    parser.add_argument(
        "--savedir", type=pathlib.Path, help="path to where the data will be saved"
    )
    parser.add_argument(
        "--gencodefile",
        type=pathlib.Path,
        help="path to where the GENCODE data in csv is saved.",
    )
    parser.add_argument(
        "--samplecol", type=str, help="Name of the column with batch (or sample) index."
    )
    parser.add_argument(
        "--celltypecol", type=str, help="Name of the column with cell types"
    )
    parser.add_argument(
        "--malignantcells",
        nargs="+",
        help="Names of the malignant cells in the celltypecol",
    )
    parser.add_argument(
        "--tmecells",
        nargs="+",
        help="Names of the tme cells in the celltypecol",
    )
    parser.add_argument(
        "--nmintme",
        type=int,
        default=50,
        help="Minimum number of TME cells needed to run infercnv",
    )
    parser.add_argument(
        "--nmintmerefcat",
        type=int,
        default=10,
        help="Minimum number of TME cells in a specific group for it to be added for infercnv",
    )
    parser.add_argument(
        "--leidenres",
        type=float,
        default=3,
        help="Resolution used for Leiden clustering in CNV space (the bigger, the more clusters)",
    )

    return parser


def first_qc_filter(adata: sc.AnnData, sample_col: str, datatype: str) -> sc.AnnData:
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True)

    def is_outlier(adata, metric: str, nmads: int):
        M = adata.obs[metric]
        outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
            np.median(M) + nmads * median_abs_deviation(M) < M
        )
        return outlier

    adatas = []
    for patient in tqdm(adata.obs[sample_col].unique()):
        patadata = adata[adata.obs[sample_col] == patient].copy()
        patadata.obs["outlier"] = (
            is_outlier(patadata, "log1p_total_counts", 5)
            | is_outlier(patadata, "log1p_n_genes_by_counts", 5)
            | is_outlier(patadata, "pct_counts_in_top_50_genes", 5)
        )

        adatas.append(patadata)

    adata = adatas[0].concatenate(*adatas[1:])
    if datatype == "10x":
        adata = adata[
            (adata.obs.total_counts > 1500)
            & (adata.obs.total_counts < 50000)
            & (adata.obs.n_genes_by_counts > 500)
            & (~adata.obs.outlier)
        ].copy()
    else:
        adata = adata[(~adata.obs.outlier)].copy()
    return adata


def remove_doublets(adata: sc.AnnData, sample_col: str) -> sc.AnnData:
    sc.external.pp.scrublet(adata, batch_key=sample_col)
    adata = adata[~adata.obs.predicted_doublet].copy()
    return adata


def infercnv_and_assign(
    patadata: sc.AnnData,
    celltype_col: str,
    tme_cells: List[str],
    malignant_cells: List[str],
    n_min_tme: int = 50,
    n_min_tme_refcat: int = 10,
    leiden_res: float = 3,
) -> Optional[sc.AnnData]:

    vc_tme = patadata.obs[celltype_col][
        patadata.obs[celltype_col].isin(tme_cells)
    ].value_counts()
    if vc_tme.sum() < n_min_tme:
        LOGGER.info(f"Patient has less than {n_min_tme} TME cells; skipping")
        return None

    tme_to_use = vc_tme[vc_tme >= n_min_tme_refcat].index.to_numpy()
    if len(tme_to_use) == 0:
        LOGGER.info(
            f"Patient has no TME groups that has more than {n_min_tme_refcat} TME cells; skipping"
        )
        return None

    LOGGER.info("Running CNV inference...")
    cnvpy.tl.infercnv(patadata, reference_key=celltype_col, reference_cat=tme_to_use)

    cnvpy.tl.pca(patadata)
    cnvpy.pp.neighbors(patadata)
    cnvpy.tl.leiden(patadata, resolution=leiden_res)

    leiden_composition = (
        patadata.obs[["cnv_leiden", celltype_col]].value_counts().unstack().fillna(0)
    )
    healthy_sum = leiden_composition[
        leiden_composition.columns.intersection(tme_cells)
    ].sum(axis=1)
    full_sum = leiden_composition[
        leiden_composition.columns.intersection(np.append(tme_cells, malignant_cells))
    ].sum(axis=1)
    cluster_assignments = (
        ((healthy_sum / full_sum) > 0.5)
        .replace({True: "TME", False: "Malignant"})
        .to_dict()
    )

    LOGGER.info("Assigning malignant status and cleaning...")
    patadata.obs["cnv_malignant_status"] = patadata.obs.cnv_leiden.replace(
        cluster_assignments
    )

    patadata.obs["cleaned_celltype"] = patadata.obs[celltype_col].astype(str).copy()
    patadata.obs.loc[
        (patadata.obs[celltype_col].isin(tme_cells))
        & (patadata.obs["cnv_malignant_status"] == "Malignant"),
        "cleaned_celltype",
    ] = "Uncertain"
    patadata.obs.loc[
        (patadata.obs[celltype_col].isin(malignant_cells))
        & (patadata.obs["cnv_malignant_status"] == "TME"),
        "cleaned_celltype",
    ] = "Uncertain"

    return patadata


def get_filtered_adata(
    adata: sc.AnnData,
    sample_col: str,
    datatype: str,
    celltype_col: str,
    tme_cells: List,
    malignant_cells: List,
    n_min_tme: int,
    n_min_tme_refcat: int,
    leiden_res: float,
) -> sc.AnnData:
    adatas = []
    for patient in tqdm(adata.obs[sample_col].unique()):
        print(patient)
        patadata = adata[adata.obs[sample_col] == patient].copy()

        if datatype == "10x":
            patadata.layers["counts"] = patadata.X.copy()
            sc.pp.normalize_total(patadata)
            sc.pp.log1p(patadata)
        else:
            patadata.X = patadata.layers["TPM"].copy()
            sc.pp.log1p(patadata)

        patadata = infercnv_and_assign(
            patadata=patadata,
            celltype_col=celltype_col,
            tme_cells=tme_cells,
            malignant_cells=malignant_cells,
            n_min_tme=n_min_tme,
            n_min_tme_refcat=n_min_tme_refcat,
            leiden_res=leiden_res,
        )
        if patadata is not None:
            adatas.append(patadata)

    filtered = adatas[0].concatenate(*adatas[1:])
    return filtered


def main() -> None:
    # Read the CLI arguments
    parser = create_parser()
    args = parser.parse_args()

    print("TME cells", args.tmecells)
    print("Malignant cells", args.malignantcells)

    LOGGER.info("Loading data...")
    adata = sc.read_h5ad(args.data)
    adata.var_names_make_unique()

    LOGGER.info("Loading gencode info...")
    gencode = pd.read_csv(args.gencodefile, index_col=0).set_index("gene_name")
    gencode = gencode.loc[~gencode.index.duplicated()].rename(
        columns={"seqname": "chromosome"}
    )

    if args.datatype == "tpm":
        LOGGER.info("Using TPM; rounding to nearest integer...")
        adata.layers["TPM"] = adata.X.copy()
        adata.X = adata.X.astype(int)
    else:
        LOGGER.info("Using 10X")
    n_cells = adata.n_obs
    LOGGER.info(f"There are n={n_cells}.")
    LOGGER.info("First QC pass...")
    adata = first_qc_filter(
        adata=adata, sample_col=args.samplecol, datatype=args.datatype
    )
    LOGGER.info(
        f"Filtered out n={n_cells - adata.n_obs} cells, n={adata.n_obs} cells left."
    )
    n_cells = adata.n_obs
    LOGGER.info("Removing doublets...")
    adata = remove_doublets(adata=adata, sample_col=args.samplecol)
    LOGGER.info(
        f"Filtered out n={n_cells - adata.n_obs} cells, n={adata.n_obs} cells left."
    )

    LOGGER.info("Adding GENCODE info...")
    adata.var = pd.concat(
        [adata.var, gencode.loc[gencode.index.intersection(adata.var.index)]], axis=1
    )

    LOGGER.info("Starting inferCNV cleaning...")
    adata = get_filtered_adata(
        adata=adata,
        sample_col=args.samplecol,
        datatype=args.datatype,
        celltype_col=args.celltypecol,
        tme_cells=args.tmecells,
        malignant_cells=args.malignantcells,
        n_min_tme=args.nmintme,
        n_min_tme_refcat=args.nmintmerefcat,
        leiden_res=args.leidenres,
    )

    os.makedirs(args.savedir, exist_ok=True)

    LOGGER.info("Saving cleaned dataset...")
    adata.write_h5ad(args.savedir / "filtered_adata.h5ad")


if __name__ == "__main__":
    main()
