#!/bin/bash

AVAILABLE_GPUS="${1:-0,1,2,3}"
chmod +x gpumaestro.sh
chmod +x gpu_task.sh

# Build the task queue
declare -a TASK_QUEUE
TASKS=("Task1" "Task2" "Task3" "Task4" "Task5" "Task6")

for task in "${TASKS[@]}"; do
    TASK_QUEUE+=("$task _GPU_")  # Simple tasks with just a name
done

echo "Starting dummy workload using available GPUs: $AVAILABLE_GPUS"

# Pass GPU list, script name, and task queue to gpumaestro.sh
./gpumaestro.sh "$AVAILABLE_GPUS" "./gpu_task.sh" "${TASK_QUEUE[@]}"