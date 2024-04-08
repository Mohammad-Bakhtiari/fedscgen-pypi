#!/bin/bash

chmod +x fedscgen_tuning.sh

echo "hyperparameter tuning for including all cell types and batches"
DATASETS=(HumanDendriticCells MouseCellAtlas HumanPancreas PBMC CellLine MouseRetina MouseHematopoieticStemProgenitorCells)
for ds in "${DATASETS[@]}";do
  echo "Running hyperparameter tuning for $ds"
  batches="0,1"
  n_clients=2
  if [ "$ds" == "HumanPancreas" ]; then
    batches="0,1,2,3,4"
    n_clients=5
  elif [ "$ds" == "CellLine" ]; then
    batches="0,1,2"
    n_clients=3
  fi
  ./fedscgen_tuning.sh "${ds}.h5ad" "$n_clients" "$batches"
done