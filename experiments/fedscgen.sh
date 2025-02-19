#!/bin/bash


# Accepting arguments
H5AD_FILE="$1"
REMOVE_CELL_TYPES="$2"
COMBINE="$3"
DROP="$4"
BATCH_OUT_VALUES=($5)
N_CLIENTS_VALUES=($6)
BATCHES=$7
GPU="${8:-1}"
SMPC="${9:-false}"
ROUND="${10:-8}"
EPOCH="${11:-2}"
BATCH_SIZE="${12:-50}"
SNAPSHOT="${13:-flase}"
TUNING="${14:-false}"


# DATASET is H5AD_FILE without the extension
DATASET=$(echo "$H5AD_FILE" | cut -f 1 -d '.')



echo "Running arguments: h5ad file: $H5AD_FILE, remove cell types: $REMOVE_CELL_TYPES, combine: $COMBINE, drop: $DROP, batch out values: ${BATCH_OUT_VALUES[@]}, n clients: ${N_CLIENTS_VALUES[@]}, batches: $BATCHES, gpu: $GPU, rounds: $ROUND, epoch: $EPOCH, batch size: $BATCH_SIZE, snapshot: $SNAPSHOT, tuning: $TUNING, smpc: $SMPC"

# Setting up other variables based on the flags
if [ "$COMBINE" = "true" ]; then
  INCLUSION="combined"
elif [ "$DROP" = "true" ]; then
  INCLUSION="dropped"
else
  INCLUSION="all"
fi

if [ "$SNAPSHOT" = "true" ]; then
  snapshot_flag="--per_round_snapshots"
fi


root_dir="$(dirname "$PWD")"
raw="${root_dir}/data/datasets/${H5AD_FILE}"
output_path="${root_dir}/results/scgen/federated"
if [ "$SMPC" = "true" ]; then
    output_path="${output_path}/smpc"
fi
output_path="${output_path}/${DATASET}/${INCLUSION}"
if [ "$TUNING" = "true" ]; then
  output="${root_dir}/results/scgen/federated/param-tuning"
  if [ "$SMPC" = "true" ]; then
    output="${output}/smpc"
  fi
  output="${output}/${DATASET}/E${EPOCH}"
fi
combine_flag=""
if [ "$COMBINE" = "true" ]; then
    combine_flag="--combine"
fi
smpc_flag=""
if [ "$SMPC" = "true" ]; then
    smpc_flag="--smpc"
fi

for i in "${!BATCH_OUT_VALUES[@]}"; do
    batch_out=${BATCH_OUT_VALUES[$i]}
    n_clients=${N_CLIENTS_VALUES[$i]}
    if [ "$TUNING" = "false" ]; then
        output="$output_path/BO${batch_out}-C${n_clients}"
    fi
    mkdir -p "${output}"
    echo "Running $batch_out batch out for $n_clients clients"
    export CUBLAS_WORKSPACE_CONFIG=:4096:8
    python3 "${root_dir}/scripts/fedscgen_.py" \
        --init_model_path "${root_dir}/models/${DATASET}" \
        --adata "$raw" \
        --output  "$output"\
        --epoch $EPOCH \
        --cell_key "cell_type" \
        --batch_key "batch" \
        --batches "$BATCHES" \
        --lr 0.001 \
        --batch_size "$BATCH_SIZE" \
        --hidden_size "800,800" \
        --z_dim 10 \
        --early_stopping_kwargs '{"early_stopping_metric": "val_loss", "patience": 20, "threshold": 0, "reduce_lr": True, "lr_patience": 13, "lr_factor": 0.1}' \
        --batch_out "$batch_out" \
        --n_clients "$n_clients" \
        --remove_cell_types "$REMOVE_CELL_TYPES" \
        --gpu "$GPU" \
        --n_rounds $ROUND   \
        $combine_flag \
        $snapshot_flag \
        $smpc_flag
    if [ "$SNAPSHOT" = "true" ]; then
        for corrected in "$output"/*.h5ad; do
          echo -e "\e[33mPCA on $corrected\e[0m \n "
          python "${root_dir}/scripts/pca_reduction_simple.py" --path "$corrected" \
            --n_components 20 \
            --output_dir "$corrected"
        done
    else
      while IFS= read -r -d '' corrected
      do
        echo -e "\e[33mPCA on $corrected\e[0m \n "
        python "${root_dir}/scripts/pca_reduction_simple.py" --path "$corrected" \
          --n_components 20 \
          --output_dir "$corrected"
      done < <(find "$output" -name '*.h5ad' -print0)
    fi
done