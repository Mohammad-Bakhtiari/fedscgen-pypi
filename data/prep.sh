#!/bin/bash

datasets=(dataset1 dataset2 dataset4 dataset5 dataset6 dataset7 dataset8 dataset10)
root_dir="$(dirname "$PWD")"

# loop through datasets
for ds in "${datasets[@]}";do
  python prep_benchmark_ds.py --root_path "${root_dir}/data/${ds}"  --dataset "$ds"
  h5ad_file=$(find "${root_dir}/data/${ds}" -name "*.h5ad")
  python "${root_dir}/scripts/pca_reduction_simple.py" --path "$h5ad_file" \
    --n_components 20 \
    --output_dir "$h5ad_file"
  mv "$h5ad_file" "datasets/"
done