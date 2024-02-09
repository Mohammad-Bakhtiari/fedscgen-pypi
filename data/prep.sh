#!/bin/bash

datasets=(dataset1 dataset2 dataset4 dataset5 dataset6 dataset7 dataset8 dataset10)

# loop through datasets
for ds in "${datasets[@]}";do
  python prep_benchmark_ds.py --root_path "/home/bba1658/FedScGen/data/${ds}"  --dataset "$ds"
done
