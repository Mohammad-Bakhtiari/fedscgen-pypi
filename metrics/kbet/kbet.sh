scenario=$1

# scenarios ==> "datasets" "batch-out"
parent_dir="$(dirname "$PWD")"
root_dir="$(dirname "$parent_dir")"

INCLUSION=("all" "dropped" "combined")
DATASETS=("HumanDendriticCells" "MouseCellAtlas" "HumanPancreas" "PBMC" "CellLine" "MouseRetina" "MouseBrain" "MouseHematopoieticStemProgenitorCells")

# Store all commands to run later
declare -a commands
missing_files=0


if [[ "${scenario}" == "approach" ]]; then
  for inclusion in "${INCLUSION[@]}"; do
    for dataset in "${DATASETS[@]}"; do
      if [[ ( "$dataset" == "HumanDendriticCells" || "$dataset" == "CellLine" ) && "$inclusion" != "all" ]]; then
        continue
      fi
      n_clients=2
      case "$dataset" in
        "HumanPancreas") n_clients=5 ;;
        "CellLine") n_clients=3 ;;
      esac

      output_dir="${root_dir}/results/fedscgen/${dataset}/${inclusion}/BO0-C${n_clients}"
      scgen="${root_dir}/results/scgen/${dataset}/${inclusion}/corrected.h5ad"
      fedscgen="${root_dir}/results/fedscgen/${dataset}/${inclusion}/BO0-C${n_clients}/fed_corrected.h5ad"

      # Check for missing files
      if [[ ! -f "$scgen" ]]; then
        echo "ERROR: Missing scGen corrected file: $scgen"
        missing_files=1
      fi
      if [[ ! -f "$fedscgen" ]]; then
        echo "ERROR: Missing FedscGen corrected file: $fedscgen"
        missing_files=1
      fi

      # Store command
      commands+=("Rscript kbet_calculator.R \"${scgen}\" \"${fedscgen}\" \"${output_dir}\" \"${dataset}\" &")
    done
  done
fi

if [[ "${scenario}" == "batch-out" ]]; then
  BATCH_OUT=(0 1 2 3 4)
  for bo in "${BATCH_OUT[@]}"; do
    scgen="${root_dir}/results/scgen/HumanPancreas/all/corrected.h5ad"
    fedscgen="${root_dir}/results/fedscgen/HumanPancreas/all/BO1-C4/${bo}/fed_corrected_with_new_studies.h5ad"
    output_dir="${root_dir}/results/fedscgen/HumanPancreas/all/BO1-C4/${bo}"

    # Check for missing files
    if [[ ! -f "$scgen" ]]; then
      echo "ERROR: Missing scGen corrected file: $scgen"
      missing_files=1
    fi
    if [[ ! -f "$fedscgen" ]]; then
      echo "ERROR: Missing FedscGen corrected file: $fedscgen"
      missing_files=1
    fi

    # Store command
    commands+=("Rscript kbet_calculator.R \"${scgen}\" \"${fedscgen}\" \"${output_dir}\" \"HumanPancreas\" &")
  done
fi

if [[ $missing_files -eq 1 ]]; then
  echo "Aborting execution due to missing files."
  exit 1
fi

echo "All required files exist. Proceeding with execution."
for cmd in "${commands[@]}"; do
  echo "Running: $cmd"
  eval "$cmd"
done
wait