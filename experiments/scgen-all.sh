#!/bin/bash

AVAILABLE_GPUS="${1:-0,1,2,3}"
source ./config.sh
declare -a TASK_QUEUE

for seed in "${SEEDS[@]}"
do
  for dataset in "${DATASETS[@]}"
  do
    if [ "$dataset" == "MouseBrain" ] && [ "$seed" -ne 42 ]; then
      continue
    fi
    task_name="${dataset}-${seed}"
    task="$task_name|${dataset}.h5ad|''|false|false|_GPU_|$seed"
    TASK_QUEUE+=("$task")
  done
done

chmod +x gpumaestro.sh
chmod +x scgen.sh
./gpumaestro.sh "$AVAILABLE_GPUS" "./scgen.sh" "${TASK_QUEUE[@]}"