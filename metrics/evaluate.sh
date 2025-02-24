#!/bin/bash

root_dir="$(dirname "$PWD")"

export IGNORE_CRYPTEN=1
chmod +x benchmark.sh
echo "Running benchmark for param tuning of FedscGen all rounds and epochs"
./benchmark.sh "tuning" false

echo "Running benchmark for batchout of scGen and FedscGen"
./benchmark.sh "batch-out" false

echo "Running benchmark for all datasets"
./benchmark.sh "datasets" false


cd "${root_dir}/metrics/lisi"
chmod +x evaluate.sh
echo "Running lisi for raw data, and corrected data by scGen and FedscGen"
./evaluate.sh "datasets"

cd ..
echo "plotting LISI results for raw, scGen, and FedscGen on all datasets"
python plot.py --scenario "lisi" --data_dir "${root_dir}/results/scgen/federated"

cd "${root_dir}/metrics/kbet"
chmod +x kbet.sh
echo "Running kbet for corrected data by scGen and FedscGen on all datsets"
./kbet.sh "datasets"

echo "Running kbet for corrected data by scGen and FedscGen on HumanPancreas batchout"
./kbet.sh "batch-out"

cd "${root_dir}/metrics"

python plot.py --scenario "tuning" \
  --data_dir "${root_dir}/results/scgen/federated/param-tuning"

python plot.py --scenario "batchout" \
  --bo_metrics_file "${root_dir}/results/scgen/federated/HumanPancreas/all/BO1-C4/batchout_metrics.csv" \
  --bo_kbet_dir "${root_dir}/results/scgen/federated/HumanPancreas/all/BO1-C4" \
  --all_ds_metrics_file "${root_dir}/results/scgen/federated/fed_cent_metrics-all.csv" \
  --hp_kbet_file "${root_dir}/results/scgen/federated/HumanPancreas/all/BO0-C5/kBET_summary_results.csv" \
  --output_dir "${root_dir}/results/scgen/federated"

python plot.py --scenario "scenarios" \
  --data_dir "${root_dir}/results/scgen/federated"

python plot.py --scenario "datasets" \
  --data_dir "${root_dir}/results/scgen/federated"

python plot.py --scenario "kbet-diff" \
  --data_dir "${root_dir}/results/scgen/federated"


echo "plotting umaps for corrected data by scGen and FedscGen on different datasets with various rounds and epochs"
python plot_umaps.py --scenario "tuning" \
  --data_dir "${root_dir}/results/scgen/centralized" \
  --raw_data_dir "${root_dir}/data/datasets" \
  --fed_data_dir "${root_dir}/results/scgen/federated/param-tuning" \
  --output_dir "${root_dir}/results/scgen/umap/tuning" \
  --round 8 \
  --epoch 2

echo "plotting umaps for corrected data by scGen and FedscGen on HumanPancreas batchout"
python plot_umaps.py --scenario "batchout" \
  --data_dir "${root_dir}/results/scgen/centralized/HumanPancreas/all" \
  --fed_data_dir "${root_dir}/results/scgen/federated/HumanPancreas/all/BO1-C4" \
  --output_dir "${root_dir}/results/scgen/umap/batch-out" \
  --n_batches 5

echo "plotting umaps for corrected data by scGen and FedscGen on all datasets"
python plot_umaps.py --scenario "datasets" \
 --data_dir "${root_dir}/results/scgen/centralized" \
 --raw_data_dir "${root_dir}/data/datasets" \
 --fed_data_dir "${root_dir}/results/scgen/federated" \
 --output_dir "${root_dir}/results/scgen/umap/datasets"

echo "plot classification accuracy for corrected data by FedscGen and scGen"
python plot.py --scenario "classification" \
  --data_dir "${root_dir}/results/scgen"