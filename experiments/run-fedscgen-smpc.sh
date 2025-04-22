#!/bin/bash
AVAILABLE_GPUS="${1:-0,1,2,3}"
source ./config.sh

declare -a TASK_QUEUE

for seed in "${SEEDS[@]}"
do
  for dataset in "${DATASETS[@]}"
  do
    if [ "$dataset" == "MouseBrain" ]; then
      continue
    fi
    n_clients="2"
    batches="0,1"
    if [ "$dataset" == "HumanPancreas" ]; then
      n_clients="5"
      batches="0,1,2,3,4"
    elif [ "$dataset" == "CellLine" ]; then
      n_clients="3"
      batches="0,1,2"
    fi
    task_name="${dataset}-Seed[${seed}]"
    task="$task_name|${dataset}.h5ad|$n_clients|$batches|_GPU_|$seed"
    TASK_QUEUE+=("$task")
  done
done
chmod +x gpumaestro.sh
chmod +x fedscgen-smpc.sh
./gpumaestro.sh "$AVAILABLE_GPUS" "./fedscgen-smpc.sh" "${TASK_QUEUE[@]}"