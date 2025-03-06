#!/bin/bash

SMPC="${1:-false}"
NUM_GPUS="${2:-3}"
GPU=0
ROUND=8
EPOCH=2
BATCH_SIZE=50
SEEDS=(42 123 456 789 101112)  # Predefined seeds

DATASETS=(HumanDendriticCells MouseCellAtlas HumanPancreas PBMC CellLine MouseRetina MouseBrain MouseHematopoieticStemProgenitorCells)

# Set root directory
root_dir="$(dirname "$PWD")"

# Set output path based on SMPC flag
output_path_prefix="${root_dir}/results/scgen/federated"
smpc_flag=""
if [ "$SMPC" = "true" ]; then
  output_path_prefix="${output_path_prefix}/smpc"
  smpc_flag="--smpc"
fi

for seed in "${SEEDS[@]}"
do
  export FedscGen_RANDOM_SEED=$seed  # Set environment variable for seed

  for dataset in "${DATASETS[@]}"
  do
    n_clients="2"
    batches="0,1"
    if [ "$dataset" == "HumanPancreas" ]; then
      n_clients="5 4"
      batches="0,1,2,3,4"
    elif [ "$dataset" == "CellLine" ]; then
      n_clients="3"
      batches="0,1,2"
    fi
    echo -e "\e[31mRunning FedscGen for $dataset with seed $seed and n_clients=$n_clients and batches=$batches on GPU $GPU\e[0m"

    # Define paths
    H5AD_FILE="$dataset.h5ad"
    raw="${root_dir}/data/datasets/${H5AD_FILE}"
    output_path="${output_path_prefix}/${dataset}/all/seed_${seed}"

    mkdir -p "${output_path}"
    export CUBLAS_WORKSPACE_CONFIG=:4096:8

    python3 "${root_dir}/scripts/fedscgen_.py" \
        --init_model_path "${root_dir}/models/${dataset}" \
        --adata "$raw" \
        --output  "$output_path" \
        --epoch $EPOCH \
        --cell_key "cell_type" \
        --batch_key "batch" \
        --batches "$batches" \
        --lr 0.001 \
        --batch_size "$BATCH_SIZE" \
        --hidden_size "800,800" \
        --z_dim 10 \
        --early_stopping_kwargs '{"early_stopping_metric": "val_loss", "patience": 20, "threshold": 0, "reduce_lr": True, "lr_patience": 13, "lr_factor": 0.1}' \
        --batch_out "0" \
        --n_clients "$n_clients" \
        --gpu "$GPU" \
        --n_rounds $ROUND   \
        "$smpc_flag" &

    GPU=$((GPU+1))
    if [ $GPU -eq "$NUM_GPUS" ]; then
      wait  # Ensure FedscGen jobs complete before reassigning GPUs
      GPU=0
    fi
  done

  wait  # Ensure all FedscGen processes for this seed complete before PCA

  for dataset in "${DATASETS[@]}"
  do
    output_path="${output_path_prefix}/${dataset}/all/seed_${seed}"
    while IFS= read -r -d '' corrected
    do
      echo -e "\e[33mRunning PCA on $corrected\e[0m \n"
      python "${root_dir}/scripts/pca_reduction_simple.py" --path "$corrected" \
        --n_components 20 \
        --output_dir "$corrected"
    done < <(find "$output_path" -name '*.h5ad' -print0)
  done
done

wait  # Final wait to ensure all processes complete before exiting
