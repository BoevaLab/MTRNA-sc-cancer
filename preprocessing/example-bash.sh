#!/bin/bash

python preprocessing-script.py \
    /your/path/here/unfiltered_data.h5ad \
    --datatype "10x" \
    --savedir /your/path/here/savedir \
    --gencodefile /your/path/here/auxiliary_data/gencode_v41_positions.csv \
    --samplecol "sample" \
    --celltypecol "cell_type" \
    --malignantcells "Malignant" \
    --tmecells "Epithelial" "T_cell" "Endothelial" "Fibroblast" "Macrophage" "Mast" "B_cell" \
    --nmintme 50 \
    --nmintmerefcat 10 \
    --leidenres 3