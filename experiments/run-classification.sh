#!/bin/bash

# Making the classification.sh executable
chmod +x classification.sh

running=0
for model in "knn" "mlp-norm"; do
  echo "Running ${model} for all datasets"
  for dataset in "HumanDendriticCells" "MouseCellAtlas" "HumanPancreas" "PBMC" "CellLine" "MouseRetina" "MouseBrain" "MouseHematopoieticStemProgenitorCells"; do
    bo=1
    if [ $dataset == "HumanPancreas" ]; then
      bo=4
    elif [ $dataset == "CellLine" ]; then
      bo=2
    fi
    running=$((running+1))
    echo "Running ${running} ==> ${model} for ${dataset} with BO: ${bo}"
    ./classification.sh "all" "${dataset}" "${bo}" "" "${model}" >> "${model}.log" &
    if [ $running -ge 3 ]; then
      wait
      running=0
    fi
  done

done