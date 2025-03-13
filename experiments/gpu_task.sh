#!/bin/bash

TASK_NAME=$1
GPU_ID=$2

echo "[INFO] Starting $TASK_NAME on GPU $GPU_ID"
sleep $((RANDOM % 5 + 2))  # Simulate workload (2-6 seconds)
echo "[INFO] Completed $TASK_NAME on GPU $GPU_ID"
