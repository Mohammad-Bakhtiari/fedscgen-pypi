#!/bin/bash

AVAILABLE_GPUS="${1:-0,1,2,3}"
SEEDS=(42 123 456 789 101112)

declare -a TASK_QUEUE
DATASETS=(HumanDendriticCells MouseCellAtlas HumanPancreas PBMC CellLine MouseRetina MouseBrain MouseHematopoieticStemProgenitorCells)
for seed in "${SEEDS[@]}"
do
  export FedscGen_RANDOM_SEED=$seed  # Set environment variable for seed
  for dataset in "${DATASETS[@]}"
  do
    task_name="${DATASETS[$index]}-${seed}"
    task="$task_name|${dataset}.h5ad|''|false|false|_GPU_|$seed"
    TASK_QUEUE+=("$task")
  done
done
wait

chmod +x gpumaestro.sh
chmod +x scgen.sh
./gpumaestro.sh "$AVAILABLE_GPUS" "./scgen.sh" "${TASK_QUEUE[@]}"