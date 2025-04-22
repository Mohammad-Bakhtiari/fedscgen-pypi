#!/bin/bash

source ./config.sh
declare -a TASK_QUEUE

for index in "${!DROPPED_DATASETS[@]}"
do
  for inclusion in dropped combined
  do
    dataset="${DROPPED_DATASETS[$index]}"
    combined=$([ "$inclusion" == "combined" ] && echo true || echo false)
    dropped=$([ "$inclusion" == "dropped" ] && echo true || echo false)
    task_name="${dataset}-${inclusion}"
    task="$task_name|${dataset}.h5ad|${DROPPED_CELLTYPES[$index]:-''}|$combined|$dropped|_GPU_"
    TASK_QUEUE+=("$task")
  done
done

chmod +x gpumaestro.sh
chmod +x scgen.sh
./gpumaestro.sh "$AVAILABLE_GPUS" "./scgen.sh" "${TASK_QUEUE[@]}"