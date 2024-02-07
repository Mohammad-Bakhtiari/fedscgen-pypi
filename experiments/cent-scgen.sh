#!/bin/bash

# Accepting arguments
DATASET="$1"
H5AD_FILE="$2"
REMOVE_CELL_TYPES="$3"
COMBINE="$4"
DROP="$5"
GPU="${6:-1}"
BATCH_SIZE="${7:-50}"
Z_DIM="${8:-10}"

# Setting up other variables based on the flags
if [ "$COMBINE" = "true" ]; then
  TARGET_FOLDER="combined"
elif [ "$DROP" = "true" ]; then
  TARGET_FOLDER="dropped"
else
  TARGET_FOLDER="all"
fi

root_dir="$(dirname "$PWD")"
raw="${root_dir}/datasets/datasets/${H5AD_FILE}"
output_path="${root_dir}/results/scgen/centralized/${DATASET}/${TARGET_FOLDER}"

combine_flag=""
if [ "$COMBINE" = "true" ]; then
    combine_flag="--combine"
fi

# Running python scripts
python "${root_dir}/centralized_scgen.py" \
--model_path "${root_dir}/models/centralized/${DATASET}" \
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