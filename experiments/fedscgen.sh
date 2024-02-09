#!/bin/bash


# Accepting arguments
H5AD_FILE="$1"
REMOVE_CELL_TYPES="$2"
COMBINE="$3"
DROP="$4"
BATCH_OUT_VALUES=($5)
N_CLIENTS_VALUES=($6)
BATCHES=$7
GPU="${8:-1}"
BATCH_SIZE="${9:-50}"


# DATASET is H5AD_FILE without the extension
DATASET=$(echo "$H5AD_FILE" | cut -f 1 -d '.')



echo $BATCHES

# Setting up other variables based on the flags
if [ "$COMBINE" = "true" ]; then
  TARGET_FOLDER="combined"
elif [ "$DROP" = "true" ]; then
  TARGET_FOLDER="dropped"
else
  TARGET_FOLDER="all"
fi

root_dir="$(dirname "$PWD")"
raw="${root_dir}/data/datasets/${H5AD_FILE}"
output_path="${root_dir}/results/scgen/federated/${DATASET}/${TARGET_FOLDER}"

combine_flag=""
if [ "$COMBINE" = "true" ]; then
    combine_flag="--combine"
fi

for i in "${!BATCH_OUT_VALUES[@]}"; do
    batch_out=${BATCH_OUT_VALUES[$i]}
    n_clients=${N_CLIENTS_VALUES[$i]}
    output="$output_path/BO${batch_out}-C${n_clients}"
    mkdir -p "${output}"
    echo "Running $batch_out batch out for $n_clients clients"

    python3 "${root_dir}/scripts/fedscgen_.py" \
        --init_model_path "${root_dir}/models/centralized/${DATASET}" \
        --adata "$raw" \
        --output  "$output"\
        --epoch 5 \
        --cell_key "cell_type" \
        --batch_key "batch" \
        --batches "$BATCHES" \
        --lr 0.001 \
        --batch_size "$BATCH_SIZE" \
        --hidden_size "800,800" \
        --z_dim 10 \
        --early_stopping_kwargs '{"early_stopping_metric": "val_loss", "patience": 20, "threshold": 0, "reduce_lr": True, "lr_patience": 13, "lr_factor": 0.1}' \
        --batch_out "$batch_out" \
        --n_clients "$n_clients" \
        --remove_cell_types "$REMOVE_CELL_TYPES" \
        --gpu "$GPU" \
        --n_rounds 10 \
        $combine_flag >> "${output}/logs"
done