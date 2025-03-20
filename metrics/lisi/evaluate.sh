#!/bin/bash

parent_dir="$(dirname "$PWD")"
root_dir="$(dirname "$parent_dir")"
INCLUSION_SCENARIOS=("all" "dropped" "combined")
DATASETS=("CellLine" "PBMC" "HumanPancreas" "MouseRetina" "MouseBrain" "MouseHematopoieticStemProgenitorCells" "HumanDendriticCells" "MouseCellAtlas")

# Array to store arguments for execution
declare -a commands
missing_files=0

# Step 1: **Collect arguments and check if all required files exist**
for inclusion in "${INCLUSION_SCENARIOS[@]}"; do
  for dataset in "${DATASETS[@]}"; do
    # Skip certain datasets based on conditions
    if [[ ( "$dataset" == "HumanDendriticCells" || "$dataset" == "CellLine" ) && "$inclusion" != "all" ]]; then
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

    # Check file existence
    if [[ ! -f "$raw" ]]; then
      echo "ERROR: Missing raw data file: $raw"
      missing_files=1
    fi
    if [[ ! -f "$scgen" ]]; then
      echo "ERROR: Missing scGen corrected file: $scgen"
      missing_files=1
    fi
    if [[ ! -f "$fedscgen" ]]; then
      echo "ERROR: Missing FedscGen corrected file: $fedscgen"
      missing_files=1
    fi

    # Store command for later execution
    commands+=("Rscript lisi.R \"${raw}\" \"${scgen}\" \"${fedscgen}\" \"${output_dir}\"")
  done
done

# If any file is missing, exit before execution
if [[ $missing_files -eq 1 ]]; then
  echo "Aborting execution due to missing files."
  exit 1
fi

# Step 2: **Execute stored commands**
echo "All required files exist. Proceeding with execution."
for cmd in "${commands[@]}"; do
  echo "Running: $cmd"
  eval "$cmd"
done
