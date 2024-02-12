#!/bin/bash

root_dir="$(dirname "$PWD")"
CENTRALIZED="${root_dir}/results/scgen/centralized/HumanPancreas/all/BO0"
FEDERATED="${root_dir}/results/scgen/federated/HumanPancreas/all/BO1-C4"
FEDERATED_C5="${root_dir}/results/scgen/federated/HumanPancreas/all/BO0-C5"
DIR="${root_dir}/results/scgen/batchout"

BATCH_OUT=(0 1 2 3 4)
for bo in "${BATCH_OUT[@]}";
do
  OUTPUT_DIR="${DIR}/${bo}"
  if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
  fi
  python "${root_dir}/scripts/pca_reduction_simple.py" --path "${CENTRALIZED}/${bo}/corrected.h5ad" \
  --n_components 20 \
  --output_dir "${OUTPUT_DIR}/ScGen.h5ad"
  python "${root_dir}/scripts/pca_reduction_simple.py" --path "${FEDERATED}/${bo}/corrected.h5ad" \
  --n_components 20 \
  --output_dir "${OUTPUT_DIR}/FedScGen.h5ad"
done

OUTPUT_DIR="${DIR}/snapshot"
if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi
python "${root_dir}/scripts/pca_reduction_simple.py" --path "${FEDERATED_C5}/corrected.h5ad" \
  --n_components 20 \
  --output_dir "${OUTPUT_DIR}/FedScGen-C5.h5ad"


rounds=(1 2 3 4 5 6 7 8 9 10)

for i in "${rounds[@]}";do
    python "${root_dir}/scripts/pca_reduction_simple.py" --path "${FEDERATED_C5}/corrected_${i}.h5ad" \
      --n_components 20 \
      --output_dir "${OUTPUT_DIR}/FedScGen-C5-${i}.h5ad"
done
