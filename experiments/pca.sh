#!/bin/bash

INCLUSION=$1
H5AD_FILE=$2
BATCH_OUT=$3
N_COMPONENTS=$4
COMMON_SPACE=$5

# DATASET is H5AD_FILE without the extension
DATASET=$(echo "$H5AD_FILE" | cut -f 1 -d '.')

root_dir="$(dirname "$PWD")"
RAW="${root_dir}/data/datasets/${H5AD_FILE}"
CENTRALIZED="${root_dir}/results/scgen/centralized/${DATASET}/${INCLUSION}/corrected.h5ad"
FEDERATED="${root_dir}/results/scgen/federated/${DATASET}/${INCLUSION}/${BATCH_OUT}/corrected.h5ad"
OUTPUT_DIR="${root_dir}/results/scgen/eval/${DATASET}/${INCLUSION}"
if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
fi
common_space_flag=""
if [ "$COMMON_SPACE" = "true" ]; then
    common_space_flag="--common_space"
fi
python "${root_dir}/scripts/pca_reduction.py" --raw "${RAW}" \
--centralized "${CENTRALIZED}" \
--federated "${FEDERATED}" \
--n_components "${N_COMPONENTS}" \
--output_dir "${OUTPUT_DIR}" \
$common_space_flag
