#!/bin/bash

AVAILABLE_GPUS="${1:-0,1,2,3}"
SMPC="${2:-"--smpc"}"

declare -a TASK_QUEUE
echo "hyperparameter tuning for including all cell types and batches"
DATASETS=(HumanDendriticCells MouseCellAtlas HumanPancreas PBMC CellLine MouseRetina MouseHematopoieticStemProgenitorCells)
for ds in "${DATASETS[@]}";do
  echo -e "\e[31mRunning hyperparameter tuning for $ds\e[0m"
  batches="0,1"
  n_clients=2
  if [ "$ds" == "HumanPancreas" ]; then
    batches="0,1,2,3,4"
    n_clients=5
  elif [ "$ds" == "CellLine" ]; then
    batches="0,1,2"
    n_clients=3
  fi
  n_rounds=10
  for epoch in {1..10}; do
        task_name="${ds}-E${epoch}"
        task="$task_name|$ds.h5ad|''|false|false|0|$n_clients|$batches|_GPU_|$n_rounds|$epoch|50|true|true"
        TASK_QUEUE+=("$task")
  done
done
script_name="fedscgen.sh"
[ "$SMPC" == "--smpc" ] && script_name="fedscgen-smpc.sh"
chmod +x gpumaestro.sh
chmod +x $script_name
./gpumaestro.sh "$AVAILABLE_GPUS" "./${script_name}" "${TASK_QUEUE[@]}"