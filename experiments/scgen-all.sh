#!/bin/bash

NUM_GPUS="${1:-3}"
GPU=0
SEEDS=(42 123 456 789 101112)  # Predefined seeds

# Making the scgen.sh executable
chmod +x scgen.sh

DATASETS=(HumanDendriticCells MouseCellAtlas HumanPancreas PBMC CellLine MouseRetina MouseBrain MouseHematopoieticStemProgenitorCells)
INCLUSION="all"
root_dir="$(dirname "$PWD")"

for seed in "${SEEDS[@]}"
do
  export FedscGen_RANDOM_SEED=$seed  # Set environment variable for seed
  for dataset in "${DATASETS[@]}"
  do
    echo -e "\e[31mRunning scGen for $dataset with seed $seed on GPU $GPU\e[0m"

    # Define paths and parameters
    raw="${root_dir}/data/datasets/${dataset}.h5ad"
    output_path="${root_dir}/results/scgen/centralized/${dataset}/${INCLUSION}/seed_${seed}"  # Include seed in output path

    mkdir -p "${output_path}"
    export CUBLAS_WORKSPACE_CONFIG=:4096:8

    python "${root_dir}/scripts/scgen.py" \
    --model_path "${root_dir}/models/${dataset}" \
    --data_path "$raw" \
    --output_path "$output_path" \
    --epoch 100 \
    --batch_key "batch" \
    --cell_key "cell_type" \
    --z_dim 10 \
    --hidden_layers_sizes "800,800" \
    --batch_size 50 \
    --early_stopping_kwargs "{'early_stopping_metric': 'val_loss', 'patience': 20, 'threshold': 0, 'reduce_lr': True, 'lr_patience': 13, 'lr_factor': 0.1}" \
    --gpu "$GPU" &

    GPU=$((GPU+1))
    if [ $GPU -eq "$NUM_GPUS" ]; then
      wait  # Ensure scGen jobs complete before reassigning GPUs
      GPU=0
    fi
  done

  wait  # Ensure all scGen processes for this seed complete before PCA reduction

  for dataset in "${DATASETS[@]}"
  do
    corrected_file="${root_dir}/results/scgen/centralized/${dataset}/${INCLUSION}/seed_${seed}/corrected.h5ad"
    echo -e "\e[32mRunning PCA reduction for $dataset with seed $seed\e[0m"
    python "${root_dir}/scripts/pca_reduction_simple.py" --path "$corrected_file" \
      --n_components 20 \
      --output_dir "$corrected_file"  # **Runs PCA one at a time**
  done

done
wait  # Final wait to ensure all processes complete before exiting