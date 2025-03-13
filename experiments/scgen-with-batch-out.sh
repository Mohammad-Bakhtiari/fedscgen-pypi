#!/bin/bash


# Accepting arguments
H5AD_FILE="$1"
REMOVE_CELL_TYPES="$2"
COMBINE="$3"
DROP="$4"
BATCH_OUT_VALUES=($5)
GPU="${6:-1}"
BATCH_SIZE="${7:-50}"
Z_DIM="${8:-10}"


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


OLD_IFS="$IFS"
IFS=',' read -r -a BATCH_OUT_VALUES <<< "$5"
IFS="$OLD_IFS"

root_dir="$(dirname "$PWD")"
raw="${root_dir}/data/datasets/${H5AD_FILE}"
output_path_ds_inclusion="${root_dir}/results/scgen/${DATASET}/${INCLUSION}"

combine_flag=""
if [ "$COMBINE" = "true" ]; then
    combine_flag="--combine"
fi
echo "Batch out values: ${BATCH_OUT_VALUES[@]}"
for batch_out in "${!BATCH_OUT_VALUES[@]}"; do
  output_path="${output_path_ds_inclusion}/BO${batch_out}"
  mkdir -p "${output_path}"
  echo "Running scgen with batch out: $batch_out"
  export CUBLAS_WORKSPACE_CONFIG=:4096:8
  python "${root_dir}/scripts/scgen.py" \
  --model_path "${root_dir}/models/centralized/${DATASET}" \
  --data_path "$raw" \
  --output_path "$output_path" \
  --epoch 100 \
  --batch_key "batch" \
  --cell_key "cell_type" \
  --z_dim "$Z_DIM" \
  --hidden_layers_sizes "800,800" \
  --batch_size "$BATCH_SIZE" \
  --target_batch "$batch_out" \
  --remove_cell_types "$REMOVE_CELL_TYPES" \
  --early_stopping_kwargs "{'early_stopping_metric': 'val_loss', 'patience': 20, 'threshold': 0, 'reduce_lr': True, 'lr_patience': 13, 'lr_factor': 0.1}" \
  --gpu "$GPU" \
  $combine_flag

  # Running PCA reduction on results# Running PCA reduction on results
  corrected_file="${output_path}/corrected.h5ad"
  python "${root_dir}/scripts/pca_reduction_simple.py" --path "$corrected_file" \
  --n_components 20 \
  --output_dir "$corrected_file"
done
