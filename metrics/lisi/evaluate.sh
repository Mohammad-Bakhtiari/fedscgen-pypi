#!/bin/bash

INCLUSION=$1

parent_dir="$(dirname "$PWD")"
root_dir="$(dirname "$parent_dir")"

DATASETS=("CellLine" "PBMC" "HumanPancreas" "MouseRetina" "MouseBrain" "MouseHematopoieticStemProgenitorCells" "HumanDendriticCells" "MouseCellAtlas")

# Array to store arguments for execution
declare -a commands
missing_files=0

for dataset in "${DATASETS[@]}"; do
  # Skip certain datasets based on conditions
  if [[ ( "$dataset" == "HumanDendriticCells" || "$dataset" == "CellLine" ) && "$INCLUSION" != "all" ]]; then
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


  raw="${root_dir}/data/datasets/${dataset}.h5ad"
  scgen="${root_dir}/results/scgen/${dataset}/${INCLUSION}/corrected.h5ad"
  fedscgen="${root_dir}/results/fedscgen/${dataset}/${INCLUSION}/BO0-C${n_clients}/fed_corrected.h5ad"
  fedscgen_smpc="${root_dir}/results/fedscgen-smpc/${dataset}/fed_corrected.h5ad"

  # Check file existence in a loop
  for file in "$raw" "$scgen" "$fedscgen"; do
    if [[ ! -f "$file" ]]; then
      echo "ERROR: Missing file: $file"
      missing_files=1
    fi
  done

  # Only check fedscgen_smpc if it's not empty

  if [[ "$dataset" == "MouseBrain" ]]; then
    fedscgen_smpc="none"
  elif [[ "$INCLUSION" == "all" && ! -f "$fedscgen_smpc" ]]; then
    echo "ERROR: Missing file: $fedscgen_smpc"
    missing_files=1
  fi

  output_dir="${root_dir}/results/fedscgen/${dataset}/${INCLUSION}/BO0-C${n_clients}"
  # Store command for later execution
  commands+=("Rscript lisi.R \"${raw}\" \"${scgen}\" \"${fedscgen}\" \"${fedscgen_smpc}\" \"${output_dir}\"")
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
