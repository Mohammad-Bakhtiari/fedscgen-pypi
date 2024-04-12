#!/bin/bash

scenario=$1

scenarios=("datasets", "batch-out")
parent_dir="$(dirname "$PWD")"
root_dir="$(dirname "$parent_dir")"

INCLUSION=("all" "dropped" "combined")
DATASETS=("HumanDendriticCells" "MouseCellAtlas" "HumanPancreas" "PBMC" "CellLine" "MouseRetina" "MouseBrain" "MouseHematopoieticStemProgenitorCells")
DROPPED_DATASETS=("MouseCellAtlas" "HumanPancreas" "PBMC" "MouseRetina" "MouseBrain" "MouseHematopoieticStemProgenitorCells")



if [[ "${scenario}" == "datasets" ]]; then
  # All datasets all scenarios
  for T in "${INCLUSION[@]}"; do
    if [[ "$T" != "all" ]]; then
      DATASETS=("${DROPPED_DATASETS[@]}")
    fi

    # Loop through each dataset
    for index in "${!DATASETS[@]}"; do
      n_clients=2
      if [[ "${DATASETS[$index]}" == "HumanPancreas" ]]; then
        n_clients=5
      elif [[ "${DATASETS[$index]}" == "CellLine" ]]; then
        n_clients=3
      fi
      output_dir="${root_dir}/results/scgen/federated/${DATASETS[$index]}/${T}/BO0-C${n_clients}"
      scgen="${root_dir}/results/scgen/centralized/${DATASETS[$index]}/${T}/corrected.h5ad"
      fedscgen="${root_dir}/results/scgen/federated/${DATASETS[$index]}/${T}/BO0-C${n_clients}/fed_corrected.h5ad"
      echo "Running kbet_calculator.R for ${DATASETS[$index]} ${T}"
      Rscript kbet_calculator.R "${scgen}" "${fedscgen}" "${output_dir}" "${DATASETS[$index]}"
    done
  done
fi



# Batch out scenarios for HumanPancreas
if [[ "${scenario}" == "batch-out" ]]; then
  BATCH_OUT=(0 1 2 3 4)
  for bo in "${BATCH_OUT[@]}";
  do
    scgen="${root_dir}/results/scgen/centralized/HumanPancreas/all/corrected.h5ad"
    fedscgen="${root_dir}/results/scgen/federated/HumanPancreas/all/BO1-C4/${bo}/fed_corrected_with_new_studies.h5ad"
    output_dir="${root_dir}/results/scgen/federated/HumanPancreas/all/BO1-C4/${bo}"
    echo "Running kbet_calculator.R for HumanPancreas batchout ${bo}"
    Rscript kbet_calculator.R "${scgen}" "${fedscgen}" "${output_dir}" "HumanPancreas"
  done
fi