#!/bin/bash

chmod +x pca.sh

INCLUSION=("all" "dropped" "combined")
DATASETS=("HumanDendriticCells.h5ad" "MouseCellAtlas.h5ad" "HumanPancreas.h5ad" "PBMC.h5ad" "CellLine.h5ad" "MouseRetina.h5ad" "MouseBrain.h5ad" "MouseHematopoieticStemProgenitorCells.h5ad")
DROPPED_DATASETS=("MouseCellAtlas.h5ad" "HumanPancreas.h5ad" "PBMC.h5ad" "MouseRetina.h5ad" "MouseBrain.h5ad" "MouseHematopoieticStemProgenitorCells.h5ad")
BATCH_OUT=("BO0-C2" "BO0-C2" "BO0-C5" "BO0-C2" "BO0-C3" "BO0-C2" "BO0-C2" "BO0-C2")
DROPPED_BATCH_OUT=("BO0-C2" "BO0-C5" "BO0-C2" "BO0-C2" "BO0-C2" "BO0-C2")
N_COMPONENTS=20
COMMON_SPACE=false
for T in "${INCLUSION[@]}"; do
  if [[ "$T" != "all" ]]; then
    DATASETS=("${DROPPED_DATASETS[@]}")
    BATCH_OUT=("${DROPPED_BATCH_OUT[@]}")
  fi
  echo "${DATASETS[@]}"
  # Loop through each dataset
  for index in "${!DATASETS[@]}"; do
    echo "PCA for $T ${DATASETS[$index]} ${BATCH_OUT[$index]}"
    ./pca.sh "$T" \
     "${DATASETS[$index]}" \
     "${BATCH_OUT[$index]}" \
     "${N_COMPONENTS}" \
     "${COMMON_SPACE}"
  done
done