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
      Rscript kbet_calculator.R "${scgen}" "${fedscgen}" "${output_dir}" "${DATASETS[$index]}" &
    done
  done
  wait
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


if [[ "${scenario}" == "smpc" ]]; then
  for dataset in "${DATASETS[@]}"; do
    declare -A scgen_files
    declare -A fedscgen_files

    while IFS= read -r file; do
      if [[ "$file" =~ seed_([0-9]+) ]]; then
        seed="${BASH_REMATCH[1]}"
        scgen_files[$seed]="$file"
      fi
    done < <(find "${root_dir}/results/scgen/centralized/${dataset}/all/" -type f -name "corrected.h5ad")

    while IFS= read -r file; do
      if [[ "$file" =~ seed_([0-9]+) ]]; then
        seed="${BASH_REMATCH[1]}"
        fedscgen_files[$seed]="$file"
      fi
    done < <(find "${root_dir}/results/scgen/federated/smpc/${dataset}/all/" -type f -name "fed_corrected.h5ad")

    # Check if every scgen file has a matching fedscgen file
    for seed in "${!scgen_files[@]}"; do
      if [[ -n "${fedscgen_files[$seed]}" ]]; then
        output_dir=$(dirname "${fedscgen_files[$seed]}")
        echo "Running kbet_calculator.R for ${dataset} with seed ${seed}"
        Rscript kbet_calculator.R "${scgen_files[$seed]}" "${fedscgen_files[$seed]}" "$output_dir" "$dataset"
      else
        echo "Warning: No matching fedscgen file found for seed ${seed} in ${dataset}" >&2
      fi
    done
  done
fi