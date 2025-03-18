#!/bin/bash


# Accepting arguments
H5AD_FILE="$1"
N_CLIENTS="$2"
BATCHES="$3"
GPU="${4:-1}"
SEED="${5:-42}"

ROUND=8
EPOCH=2


# DATASET is H5AD_FILE without the extension
DATASET=$(echo "$H5AD_FILE" | cut -f 1 -d '.')

echo "Args: h5ad file: $H5AD_FILE, n clients: ${N_CLIENTS}, batches: $BATCHES, gpu: $GPU, seed: $SEED"


root_dir="$(dirname "$PWD")"
raw="${root_dir}/data/datasets/${H5AD_FILE}"
init_model_path="${root_dir}/models/${DATASET}"
output="${root_dir}/results/fedscgen-smpc/${DATASET}"
if [ "${SEED}" -ne "42" ]; then
  output="${output}/seed_${SEED}"
  init_model_path="${root_dir}/models/${DATASET}_${SEED}"
fi

mkdir -p "${output}"
export CUBLAS_WORKSPACE_CONFIG=:4096:8
python3 "${root_dir}/scripts/fedscgen_.py" \
        --init_model_path "$init_model_path" \
        --adata "$raw" \
        --output  "$output" \
        --epoch $EPOCH \
        --cell_key "cell_type" \
        --batch_key "batch" \
        --batches "$BATCHES" \
        --lr 0.001 \
        --batch_size 50 \
        --hidden_size "800,800" \
        --z_dim 10 \
        --early_stopping_kwargs '{"early_stopping_metric": "val_loss", "patience": 20, "threshold": 0, "reduce_lr": True, "lr_patience": 13, "lr_factor": 0.1}' \
        --batch_out "0" \
        --n_clients "$N_CLIENTS" \
        --gpu "$GPU" \
        --n_rounds $ROUND   \
        --aggregation "fedavg" \
        --seed "$SEED" \
        --smpc &

while IFS= read -r -d '' corrected
do
  echo -e "\e[33mPCA on $corrected\e[0m \n "
  python "${root_dir}/scripts/pca_reduction_simple.py" --path "$corrected" \
    --n_components 20 \
    --output_dir "$corrected"
done < <(find "$output" -name '*.h5ad' -print0)