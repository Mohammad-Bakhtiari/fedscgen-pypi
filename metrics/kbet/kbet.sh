#!/bin/bash


parent_dir="$(dirname "$PWD")"
root_dir="$(dirname "$parent_dir")"

TARGET_FOLDER=("all" "dropped" "combined")
DATASETS_NUMBERS=(1 2 4 5 6 7 8 10)
DROPPED_DATASETS_NUMBERS=(2 4 5 7 8 10)

# All datasets all scenarios
for T in "${TARGET_FOLDER[@]}"; do
  if [[ "$T" != "all" ]]; then
    DATASETS_NUMBERS=("${DROPPED_DATASETS_NUMBERS[@]}")
  fi

  # Loop through each dataset
  for index in "${!DATASETS_NUMBERS[@]}"; do
    DIR_PATH="${root_dir}/results/scgen/eval/dataset${DATASETS_NUMBERS[$index]}/${T}"
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