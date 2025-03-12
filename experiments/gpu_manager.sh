#!/bin/bash

GPU_LIST="${1:-0,1,2,3}"  # Default GPUs
IFS=',' read -r -a GPUS <<< "$GPU_LIST"  # Convert to array
NUM_GPUS=${#GPUS[@]}  # Get the number of GPUs
GPU_INDEX=0  # Counter for GPU selection
export FEDSCGEN_NEXT_GPU  # Global variable to hold the next GPU
RUNNING_TASKS=0  # Counter for active tasks

# Function to set the next available GPU
get_next_gpu() {
    FEDSCGEN_NEXT_GPU=${GPUS[$GPU_INDEX]}  # Set the global FEDSCGEN_NEXT_GPU variable
    GPU_INDEX=$(( (GPU_INDEX + 1) % NUM_GPUS ))  # Increment and reset if at end
    ((RUNNING_TASKS++))  # Increment running tasks
}

# Function to wait only when all GPUs are in use
wait_for_free_gpu() {
    # Wait if the number of running tasks equals or exceeds NUM_GPUS
    if [ "$RUNNING_TASKS" -ge "$NUM_GPUS" ]; then
        wait -n  # Wait for any child process to finish
        ((RUNNING_TASKS--))  # Decrement running tasks when one finishes
    fi
}