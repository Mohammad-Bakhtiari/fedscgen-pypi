#!/bin/bash

NUM_GPUS=3
GPU=0
chmod +x fedscgen.sh

echo "hyperparameter tuning for including all cell types and batches"
DATASETS=(HumanDendriticCells MouseCellAtlas HumanPancreas PBMC CellLine MouseRetina MouseHematopoieticStemProgenitorCells)
for ds in "${DATASETS[@]}";do
  echo "\e[31mRunning hyperparameter tuning for $ds\e[0m"
  batches="0,1"
  n_clients=2
  if [ "$ds" == "HumanPancreas" ]; then
    batches="0,1,2,3,4"
    n_clients=5
  elif [ "$ds" == "CellLine" ]; then
    batches="0,1,2"
    n_clients=3
  fi
  n_rounds=10
  epoch=1
  while true; do
      ./fedscgen.sh "$ds.h5ad" "" false false "0" "$n_clients" "$batches" "$GPU" "$n_rounds" "$epoch" 50 true true &
      GPU=$((GPU+1))
      if [ $GPU -eq $NUM_GPUS ]; then
          wait
          GPU=0
      fi
      epoch=$((epoch+1))
      if [ $epoch -gt 10 ]; then
          break
      fi
      break
  done
done