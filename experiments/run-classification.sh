#!/bin/bash

N_GPUS="${1:-3}"
BATCH_SIZE=128
running=0
root_dir="$(dirname "$PWD")"

# Define datasets with associated client numbers
declare -A DATASETS=(
  ["HumanDendriticCells"]=2
  ["MouseCellAtlas"]=2
  ["HumanPancreas"]=5
  ["PBMC"]=2
  ["CellLine"]=3
  ["MouseRetina"]=2
  ["MouseBrain"]=2
  ["MouseHematopoieticStemProgenitorCells"]=2
)

for model in "knn" "mlp-norm"; do
  echo "Running ${model} for all datasets"

  for dataset in "${!DATASETS[@]}"; do
    n_clients="${DATASETS[$dataset]}"

    export CUBLAS_WORKSPACE_CONFIG=:4096:8

    for config_key in "scgen" "fedscgen" "fedscgen-smpc"; do

      if [ "$config_key" == "scgen" ]; then
        OUTPUT="${root_dir}/results/scgen/centralized/${dataset}/all/classification/${model}"
        ADATA="${root_dir}/results/scgen/centralized/${dataset}/all/corrected.h5ad"
      else
        OUTPUT="${root_dir}/results/scgen/federated"
        [ "$config_key" == "fedscgen-smpc" ] && OUTPUT="${OUTPUT}/smpc"
        ADATA="${OUTPUT}/${dataset}/all/BO0-C${n_clients}/fed_corrected.h5ad"
        OUTPUT="${OUTPUT}/${dataset}/all/BO0-C${n_clients}/classification/${model}"
      fi
      echo "Running ${config_key} with ${model} on ${dataset} (GPU:${running})"
      python "${root_dir}/scripts/classification.py" \
        --adata "$ADATA" \
        --output "$OUTPUT" \
        --cell_key "cell_type" \
        --batch_key "batch" \
        --epoch 50 \
        --lr 0.001 \
        --batch_size $BATCH_SIZE \
        --init_model_path "${root_dir}/models/classification/${dataset}/model.pth" \
        --norm_method "min_max" \
        --model "$model" \
        --hidden_size "800,200" \
        --batch_out "1" \
        --gpu "$running" &
      if [ "$model" == "knn" ]; then
        wait
      else
        running=$((running+1))
        [ $running -ge $N_GPUS ] && wait && running=0
      fi
    done
  done
done
wait