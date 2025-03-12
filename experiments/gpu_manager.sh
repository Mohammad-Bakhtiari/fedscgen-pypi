#!/bin/bash

GPU_LIST="${1:-0,1,2,3}"  # Default GPUs
IFS=',' read -r -a GPUS <<< "$GPU_LIST"  # Convert to array
NUM_GPUS=${#GPUS[@]}  # Get the number of GPUs
export GPU_INDEX=0  # Export GPU_INDEX to persist in parent shell

# Function to get the next available GPU
get_next_gpu() {
    echo "DEBUG: Current GPU_INDEX before update: $GPU_INDEX" >&2  # Debugging

    gpu=${GPUS[$GPU_INDEX]}  # Select GPU

    # Update GPU_INDEX globally
    GPU_INDEX=$(( (GPU_INDEX + 1) % NUM_GPUS ))

    echo "DEBUG: Next GPU_INDEX after update: $GPU_INDEX" >&2  # Debugging
    echo "$gpu"  # Return the selected GPU
}

# Function to wait only for this script’s processes
wait_for_free_gpu() {
    local script_pid=$$  # Get the current script's process ID
    while [ "$(pgrep -P "$script_pid" | wc -l)" -ge "$NUM_GPUS" ]; do
        wait -n  # Wait for any process started by this script
    done
}