#!/bin/bash

GPU_LIST="${1:-0,1,2,3}"  # Default GPUs
IFS=',' read -r -a GPUS <<< "$GPU_LIST"  # Convert to array
NUM_GPUS=${#GPUS[@]}  # Get the number of GPUs
GPU_INDEX=0  # Counter for GPU selection
export FEDSCGEN_NEXT_GPU  # Global variable to hold the next GPU
RUNNING_TASKS=0  # Counter for active tasks
SCRIPT_TO_RUN="$2"  # Script to execute (e.g., fedscgen.sh)
shift 2  # Remove GPU_LIST and SCRIPT_TO_RUN from args
TASK_QUEUE=("$@")  # Remaining args are the task queue

# Function to set the next available GPU
get_next_gpu() {
    FEDSCGEN_NEXT_GPU=${GPUS[$GPU_INDEX]}
    GPU_INDEX=$(( (GPU_INDEX + 1) % NUM_GPUS ))
    ((RUNNING_TASKS++))  # Increment when a task is launched
}

# Function to maintain NUM_GPUS tasks running
wait_for_free_gpu() {
    if [ "$RUNNING_TASKS" -lt "$NUM_GPUS" ]; then
        echo -e "\e[32m[GPUmaestro]: No wait needed (running: $RUNNING_TASKS, max: $NUM_GPUS)\e[0m" >&2
        return
    fi
    echo -e "\e[32m[GPUmaestro]: Waiting for a task to finish (running: $RUNNING_TASKS, max: $NUM_GPUS)\e[0m" >&2
    wait -n  # Wait for any child process to finish
    ((RUNNING_TASKS--))  # Decrement when a task finishes
}

# Main logic to keep NUM_GPUS tasks running
echo "Starting GPUmaestro with script: $SCRIPT_TO_RUN"

# Launch initial tasks to fill all GPUs
while [ "$RUNNING_TASKS" -lt "$NUM_GPUS" ] && [ ${#TASK_QUEUE[@]} -gt 0 ]; do
    task="${TASK_QUEUE[0]}"
    unset 'TASK_QUEUE[0]'
    TASK_QUEUE=("${TASK_QUEUE[@]}")  # Shift array

    # Split task string into an array using '|' as the delimiter
    IFS='|' read -r -a task_args <<< "$task"
    task_name="${task_args[0]}"
    args=("${task_args[@]:1}")

    # Assign the next GPU *before* replacing the placeholder
    get_next_gpu

    # Replace _GPU_ placeholder inside args array with the current GPU
    for i in "${!args[@]}"; do
        args[$i]="${args[$i]//_GPU_/$FEDSCGEN_NEXT_GPU}"
    done

    echo -e "\e[32m[GPUmaestro]: Running $SCRIPT_TO_RUN on task: $task_name on GPU:$FEDSCGEN_NEXT_GPU\e[0m" >&2
    "$SCRIPT_TO_RUN" "${args[@]}" &
done

# Keep NUM_GPUS tasks running by replacing finished ones
while [ ${#TASK_QUEUE[@]} -gt 0 ]; do
    wait_for_free_gpu

    task="${TASK_QUEUE[0]}"
    unset 'TASK_QUEUE[0]'
    TASK_QUEUE=("${TASK_QUEUE[@]}")  # Shift array

    # Split task string into an array using '|' as the delimiter
    IFS='|' read -r -a task_args <<< "$task"
    task_name="${task_args[0]}"
    args=("${task_args[@]:1}")

    # Assign the next GPU *before* replacing the placeholder
    get_next_gpu

    # Replace _GPU_ placeholder inside args array with the current GPU
    for i in "${!args[@]}"; do
        args[$i]="${args[$i]//_GPU_/$FEDSCGEN_NEXT_GPU}"
    done

    echo -e "\e[32m[GPUmaestro]: Running $SCRIPT_TO_RUN on task: $task_name on GPU:$FEDSCGEN_NEXT_GPU\e[0m" >&2
    "$SCRIPT_TO_RUN" "${args[@]}" &
done

# Wait for all remaining tasks to complete
wait
echo "All tasks completed!"