#!/bin/bash

AVAILABLE_GPUS="${1:-0,1,2,3}"

declare -a TASK_QUEUE
echo "hyperparameter tuning for including all cell types and batches"
DATASETS=(HumanDendriticCells MouseCellAtlas HumanPancreas PBMC CellLine MouseRetina MouseHematopoieticStemProgenitorCells)
for ds in "${DATASETS[@]}";do
  n_clients=$([ "${dataset}" == "HumanPancreas" ] && echo "5" || echo "2")
  batches=$([ "${dataset}" == "HumanPancreas" ] && echo "0,1,2,3,4" || echo "0,1")
  if [ "${dataset}" == "CellLine" ]; then
    n_clients="3"
    batches="0,1,2"
  fi
  for epoch in {1..10}; do
        task_name="${ds}-E${epoch}"
        task="$task_name|$ds.h5ad|''|false|false|0|$n_clients|$batches|_GPU_|10|$epoch|true|true"
        TASK_QUEUE+=("$task")
  done
done

chmod +x gpumaestro.sh
chmod +x fedscgen.sh
./gpumaestro.sh "$AVAILABLE_GPUS" "./fedscgen.sh" "${TASK_QUEUE[@]}"