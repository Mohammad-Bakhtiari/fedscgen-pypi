#!/bin/bash

parent_dir="$(dirname "$PWD")"
root_dir="$(dirname "$parent_dir")"
INCLUSION_SCENARIOS=("all" "dropped" "combined")
DATASETS=("CellLine" "PBMC" "HumanPancreas" "MouseRetina" "MouseBrain" "MouseHematopoieticStemProgenitorCells" "HumanDendriticCells" "MouseCellAtlas")

for inclusion in "${INCLUSION_SCENARIOS[@]}"; do
  for dataset in "${DATASETS[@]}"; do
    if [[ ("$dataset" == "HumanDendriticCells" ||  "$dataset" == "CellLine") && "$inclusion" != "all" ]]; then
      continue
    fi
    n_clients=2
    case "$dataset" in
      "HumanPancreas")
        n_clients=5
        ;;
      "CellLine")
        n_clients=3
        ;;
    esac

    output_dir="${root_dir}/results/fedscgen/${dataset}/${inclusion}/BO0-C${n_clients}"
    raw="${root_dir}/data/datasets/${dataset}.h5ad"
    scgen="${root_dir}/results/scgen/${dataset}/${inclusion}/corrected.h5ad"
    fedscgen="${root_dir}/results/fedscgen/${dataset}/${inclusion}/BO0-C${n_clients}/fed_corrected.h5ad"

    echo "Running lisi.R for ${dataset} in ${inclusion} scenario"
    Rscript lisi.R "${raw}" "${scgen}" "${fedscgen}" "${output_dir}"
  done
done
