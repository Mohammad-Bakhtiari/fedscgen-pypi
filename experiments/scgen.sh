#!/bin/bash

# Accepting arguments
H5AD_FILE="$1"
REMOVE_CELL_TYPES="$2"
COMBINE="$3"
DROP="$4"
GPU="${5:-0}"
BATCH_SIZE="${6:-50}"
Z_DIM="${7:-10}"

# DATASET is H5AD_FILE without the extension
DATASET=$(echo "$H5AD_FILE" | cut -f 1 -d '.')

# Setting up other variables based on the flags
if [ "$COMBINE" = "true" ]; then
  INCLUSION="combined"
elif [ "$DROP" = "true" ]; then
  INCLUSION="dropped"
else
  INCLUSION="all"
fi

root_dir="$(dirname "$PWD")"
raw="${root_dir}/data/datasets/${H5AD_FILE}"
output_path="${root_dir}/results/scgen/centralized/${DATASET}/${INCLUSION}"

combine_flag=""
if [ "$COMBINE" = "true" ]; then
    combine_flag="--combine"
fi

mkdir -p "${output_path}"

python "${root_dir}/scripts/scgen.py" \
--model_path "${root_dir}/models/${DATASET}" \
--data_path "$raw" \
--output_path "$output_path" \
--epoch 100 \
--batch_key "batch" \
--cell_key "cell_type" \
--z_dim "$Z_DIM" \
--hidden_layers_sizes "800,800" \
--batch_size "$BATCH_SIZE" \
--remove_cell_types "$REMOVE_CELL_TYPES" \
--early_stopping_kwargs "{'early_stopping_metric': 'val_loss', 'patience': 20, 'threshold': 0, 'reduce_lr': True, 'lr_patience': 13, 'lr_factor': 0.1}" \
--gpu "$GPU" \
$combine_flag

corrected_file="${output_path}/corrected.h5ad"
# Running PCA reduction on results
python "${root_dir}/scripts/pca_reduction_simple.py" --path "$corrected_file" \
  --n_components 20 \
  --output_dir "$corrected_file"