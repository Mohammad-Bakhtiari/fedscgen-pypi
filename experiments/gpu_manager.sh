#!/bin/bash

GPU_LIST="${1:-0,1,2,3}"  # Default GPUs
IFS=',' read -r -a GPUS <<< "$GPU_LIST"  # Convert to array
NUM_GPUS=${#GPUS[@]}  # Get the number of GPUs
GPU_INDEX=0  # Track which GPU to use
export NEXT_GPU  # Global variable to hold the next GPU

# Function to set the next available GPU
get_next_gpu() {
    echo "DEBUG: Current GPU_INDEX before update: $GPU_INDEX" >&2  # Debugging

    NEXT_GPU=${GPUS[$GPU_INDEX]}  # Set the global NEXT_GPU variable

    # Update GPU_INDEX globally
    GPU_INDEX=$(( (GPU_INDEX + 1) % NUM_GPUS ))

    echo "DEBUG: Next GPU_INDEX after update: $GPU_INDEX" >&2  # Debugging
}

# Function to wait only for this script’s processes
wait_for_free_gpu() {
    local script_pid=$$  # Get the current script's process ID
    while [ "$(pgrep -P "$script_pid" | wc -l)" -ge "$NUM_GPUS" ]; do
        wait -n  # Wait for any process started by this script
    done
}