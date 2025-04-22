#!/bin/bash

scenario=$1
plot_only=${2:-false}

plot_only_flag=""
if [ "$plot_only" = "true" ]; then
  plot_only_flag="--plot_only"
fi

# scenarios ==> "approach" "tuning"

root_dir="$(dirname "$PWD")"


# Approach
if [ "$scenario" = "approach" ]; then
  for approach in scgen fedscgen fedscgen-smpc; do
    python benchmark.py --data_dir "${root_dir}/results" \
    --approach "$approach" \
    --scenarios "approach" &
  done
  wait
fi

# Hyperparameter tuning rounds & epochs
if [ "$scenario" = "tuning" ]; then
  python benchmark.py --data_dir "/home/bba1658/FedscGen/results/fedscgen/param-tuning" --scenarios "tuning"
fi