#!/bin/bash

AVAILABLE_GPUS="${1:-0,1,2,3}"

declare -a TASK_QUEUE

DATASETS=(HumanDendriticCells MouseCellAtlas HumanPancreas PBMC CellLine MouseRetina MouseHematopoieticStemProgenitorCells MouseBrain)

for dataset in "${DATASETS[@]}"; do
  case "$dataset" in
    "HumanPancreas") n_clients="5" ;;
    "CellLine") n_clients="3" ;;
    *) n_clients="2" ;;
  esac
  for approach in "scgen" "fedscgen" "fedscgen-smpc"; do
    task_name="${dataset}-${approach}"
    task="$task_name|$approach|$dataset|_GPU_|$n_clients"
    TASK_QUEUE+=("$task")
  done
done
chmod +x gpumaestro.sh
chmod +x classification.sh
./gpumaestro.sh "$AVAILABLE_GPUS" "./classification.sh" "${TASK_QUEUE[@]}"