#!/bin/bash

GPU_LIST="${1:-0,1,2}"  # Default: GPUs 0,1,2 if not provided
IFS=',' read -r -a GPUS <<< "$GPU_LIST"
NUM_GPUS=${#GPUS[@]}
GPU_INDEX=0

# Function to get the next available GPU
get_next_gpu() {
    local gpu=${GPUS[$GPU_INDEX]}
    GPU_INDEX=$(( (GPU_INDEX + 1) % NUM_GPUS ))
    echo "$gpu"
}

# Function to wait only for this script’s processes
wait_for_free_gpu() {
    local script_pid=$$  # Get the current script's process ID
    while [ "$(pgrep -P "$script_pid" | wc -l)" -ge "$NUM_GPUS" ]; do
        wait -n  # Wait only for any process started by this script
    done
}
