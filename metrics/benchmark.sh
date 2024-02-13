#!/bin/bash


root_dir="$(dirname "$PWD")"

# datasets
for inclusion in "all" "dropped" "combined"; do
  python benchmark.py --data_dir "${root_dir}/results/scgen/eval" \
  --scenarios "datasets" \
  --inclusion "${inclusion}" \
  --cell_key "cell_type" \
  --batch_key "batch" \
  --n_components 20
done

# snapshot
#python benchmark.py --data_dir "${root_dir}/results/scgen/batchout/snapshot" \
#--scenarios "snapshots" \
#--inclusion "all" \
#--cell_key "cell_type" \
#--batch_key "batch" \
#--n_components 20 \
#--n_rounds 10

# batch_out
#python benchmark.py --data_dir "${root_dir}/results/scgen/batchout" \
#--scenarios "batch-out" \
#--inclusion all \
#--cell_key "cell_type" \
#--batch_key "batch" \
#--n_batches 5 \
#--n_components 20
