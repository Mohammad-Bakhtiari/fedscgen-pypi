#!/bin/bash

AVAILABLE_GPUS="${1:-0,1,2,3}"
source ./gpu_manager.sh "$AVAILABLE_GPUS"

chmod +x gpu_task.sh  # Ensure the GPU task script is executable

echo "Starting dummy workload using available GPUs: $AVAILABLE_GPUS"

TASKS=("Task1" "Task2" "Task3" "Task4" "Task5" "Task6")

for task in "${TASKS[@]}"; do
    echo -e "\e[31mRunning $task\e[0m"

    GPU=$(get_next_gpu)  # Get an available GPU
    echo "DEBUG: Running $task on GPU $GPU" >&2  # Debug message

    ./gpu_task.sh "$task" "$GPU" &  # Run task in background

    wait_for_free_gpu  # Ensure we don't overload GPUs
done

wait  # Ensure all tasks finish before exiting
echo "All tasks completed!"
