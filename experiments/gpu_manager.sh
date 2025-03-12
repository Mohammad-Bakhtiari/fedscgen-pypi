#!/bin/bash

GPU_LIST="${1:-0,1,2,3}"  # Default GPUs
IFS=',' read -r -a GPUS <<< "$GPU_LIST"  # Convert string to array
NUM_GPUS=${#GPUS[@]}  # Get number of GPUs
GPU_INDEX=0  # Track which GPU is next

# Function to get the next available GPU
get_next_gpu() {
    local gpu=${GPUS[$GPU_INDEX]}  # Select current GPU
    echo "DEBUG: Selected GPU: $gpu (Index: $GPU_INDEX / $NUM_GPUS)" >&2  # Debugging output

    GPU_INDEX=$(( (GPU_INDEX + 1) % NUM_GPUS ))  # Cycle index correctly
    echo "$gpu"  # Return the selected GPU
}

# Function to wait only for this script’s processes
wait_for_free_gpu() {
    local script_pid=$$  # Get the current script's process ID
    while [ "$(pgrep -P "$script_pid" | wc -l)" -ge "$NUM_GPUS" ]; do
        wait -n  # Wait for any process started by this script
    done
}
