#!/bin/bash


# Accepting arguments
H5AD_FILE="$1"
N_CLIENTS=$2
BATCHES=$3
GPU="${4:-2}"

# DATASET is H5AD_FILE without the extension
DATASET=$(echo "$H5AD_FILE" | cut -f 1 -d '.')
root_dir="$(dirname "$PWD")"
raw="${root_dir}/data/datasets/${H5AD_FILE}"
output_path="${root_dir}/results/scgen/federated/param-tuning/${DATASET}"

# Hyperparameter tuning is done only for batchout ZERO and all Inclusion and snapshot is set to true
n_rounds=10
epoch=1
while true; do
    output="$output_path/E${epoch}"
    echo -e "\e[32m Running $DATASET with epoch=$epoch \e[0m \n"
    mkdir -p "${output}"
    python3 "${root_dir}/scripts/fedscgen_.py" \
        --init_model_path "${root_dir}/models/${DATASET}" \
        --adata "$raw" \
        --output  "$output"\
        --epoch "$epoch" \
        --cell_key "cell_type" \
        --batch_key "batch" \
        --batches "$BATCHES" \
        --lr 0.001 \
        --batch_size 50 \
        --hidden_size "800,800" \
        --z_dim 10 \
        --early_stopping_kwargs '{"early_stopping_metric": "val_loss", "patience": 20, "threshold": 0, "reduce_lr": True, "lr_patience": 13, "lr_factor": 0.1}' \
        --batch_out 0 \
        --n_clients "$N_CLIENTS" \
        --remove_cell_types "" \
        --gpu "$GPU" \
        --n_rounds "$n_rounds" \
        --per_round_snapshots
    # PCA reduction for all corrected files
    for corrected in "$output"/*.h5ad; do
      echo -e "\e[33mPCA on $corrected\e[0m \n "
      python "${root_dir}/scripts/pca_reduction_simple.py" --path "$corrected" \
        --n_components 20 \
        --output_dir "$corrected"
    done
    epoch=$((epoch+1))
    if [ $epoch -gt 10 ]; then
        break
    fi
done