#!/bin/bash

CELL_TYPE_INCLUSION=$1
DATASET=$2
BATCH_OUT=$3
REMOVE_CELL_TYPES=$4
MODEL=$5

IFS=',' read -ra BATCH_OUT_OPTIONS <<< "$BATCH_OUT_OPTIONS_STRING"
#IFS=',' read -ra MODELS <<< "$MODELS"
echo $MODEL

# Setting up other variables based on the flags
root_dir="$(dirname "$PWD")"
n_clients=$((BATCH_OUT+1))

echo "Running ${DATASET} for ${n_clients} clients"



#CORRECTED="${root_dir}/results/scgen/centralized/${DATASET}/${TARGET_FOLDER}/corrected.h5ad"
CELL_KEY="cell_type"
CELL_KEY_ALL="cell_type"
BATCH_KEY="batch"
EPOCH=50
LR=0.001
BATCH_SIZE=32
NORM_METHOD="min_max"


declare -A HIDDEN_SIZES
HIDDEN_SIZE="800,200"
MODEL_REPO="${root_dir}/models/classification/${DATASET}"

if [ ! -d "$MODEL_REPO" ]; then
    mkdir -p "$MODEL_REPO"
fi
MODEL_REPO="${MODEL_REPO}/model.pth"

# Adjust output path based on the MODEL and HIDDEN_SIZE



# all-classes configurations
declare -A CONFIGS

CONFIGS[scgen]="$scgen"
CONFIGS[fedscgen]="$fedscgen|--cell_key|$CELL_KEY_ALL"

for config_key in "scgen" "fedscgen"; do
  echo "$config_key"
  if [ "$config_key" == "scgen" ]; then
    OUTPUT="${root_dir}/results/scgen/centralized/${DATASET}/all/classification/${MODEL}"
    ADATA="${root_dir}/results/scgen/centralized/${DATASET}/all/corrected.h5ad"
  elif [ "$config_key" == "fedscgen" ]; then
    OUTPUT="${root_dir}/results/scgen/federated/${DATASET}/all/BO0-C${n_clients}/classification/${MODEL}"
    ADATA="${root_dir}/results/scgen/federated/${DATASET}/all/BO0-C${n_clients}/fed_corrected.h5ad"
  fi
  echo $OUTPUT
  echo $ADATA


  python "${root_dir}/scripts/classification.py" \
      --adata "$ADATA" \
      --output "$OUTPUT" \
      --cell_key $CELL_KEY \
      --batch_key $BATCH_KEY \
      --epoch $EPOCH \
      --lr $LR \
      --batch_size $BATCH_SIZE \
      --init_model_path "$MODEL_REPO" \
      --norm_method $NORM_METHOD \
      --model "$MODEL" \
      --hidden_size "$HIDDEN_SIZE" \
      --batch_out "$BATCH_OUT"
done