#!/bin/bash

scenario=$1

scenarios=("datasets")
parent_dir="$(dirname "$PWD")"
root_dir="$(dirname "$parent_dir")"
INCLUSION_SCENARIOS=("all" "dropped" "combined")
DATASETS=("CellLine" "PBMC" "HumanPancreas" "MouseRetina" "MouseBrain" "MouseHematopoieticStemProgenitorCells" "HumanDendriticCells" "MouseCellAtlas")
DROPPED_DATASETS=("PBMC" "HumanPancreas" "MouseRetina" "MouseBrain" "MouseHematopoieticStemProgenitorCells" "MouseCellAtlas")


if [[ "$scenario" == "datasets" ]]; then
  for T in "${INCLUSION_SCENARIOS[@]}"; do
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
      raw="${root_dir}/data/datasets/${DATASETS[$index]}.h5ad"
      scgen="${root_dir}/results/scgen/centralized/${DATASETS[$index]}/${T}/corrected.h5ad"
      fedscgen="${root_dir}/results/scgen/federated/${DATASETS[$index]}/${T}/BO0-C${n_clients}/fed_corrected.h5ad"
      echo "Running lisi.R for ${DATASETS[$index]} in ${T} scenario"
      Rscript lisi.R "${raw}" "${scgen}" "${fedscgen}" "${output_dir}"
    done
  done
fi

