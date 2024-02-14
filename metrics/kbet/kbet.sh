#!/bin/bash


parent_dir="$(dirname "$PWD")"
root_dir="$(dirname "$parent_dir")"

INCLUSION=("all" "dropped" "combined")
DATASETS=("HumanDendriticCells" "MouseCellAtlas" "HumanPancreas" "PBMC" "CellLine" "MouseRetina" "MouseBrain" "MouseHematopoieticStemProgenitorCells")
DROPPED_DATASETS=("MouseCellAtlas" "HumanPancreas" "PBMC" "MouseRetina" "MouseBrain" "MouseHematopoieticStemProgenitorCells")

# All datasets all scenarios
for T in "${INCLUSION[@]}"; do
  if [[ "$T" != "all" ]]; then
    DATASETS=("${DROPPED_DATASETS[@]}")
  fi

  # Loop through each dataset
  for index in "${!DATASETS[@]}"; do
    DIR_PATH="${root_dir}/results/scgen/eval/${DATASETS[$index]}/${T}"
    echo "Running kbet_calculator.R for ${DIR_PATH}"
    Rscript kbet_calculator.R "${DIR_PATH}"
  done
done

# Batch out scenarios for HumanPancreas
BATCH_OUT=(0 1 2 3 4)
for bo in "${BATCH_OUT[@]}";
do
  DIR_PATH="${root_dir}/results/scgen/batchout/${bo}"
  echo "Running kbet_calculator.R for ${DIR_PATH}"
  Rscript kbet_calculator.R "${DIR_PATH}"
done

# Snapshots per round for zero left-out batch scenario for HumanPancreas dataset
DIR_PATH="${root_dir}/results/scgen/batchout/snapshot"
echo "Running kbet_calculator.R for ${DIR_PATH}"
Rscript kbet_calculator.R "${DIR_PATH}"