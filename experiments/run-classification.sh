#!/bin/bash

# Making the classification.sh executable
chmod +x classification.sh


BATCH_SIZE=128
running=0
N_GPUS=3
for model in "knn" "mlp-norm"; do
  echo "Running ${model} for all datasets"
  for dataset in "HumanDendriticCells" "MouseCellAtlas" "HumanPancreas" "PBMC" "CellLine" "MouseRetina" "MouseBrain" "MouseHematopoieticStemProgenitorCells"; do
    bo=1
    if [ $dataset == "HumanPancreas" ]; then
      bo=4
    elif [ $dataset == "CellLine" ]; then
      bo=2
    fi
    ./classification.sh "all" "${dataset}" "${bo}" "" "${model}" $BATCH_SIZE $running >> "${model}.log" &
    running=$((running+1))
    echo "Running ${running} ==> ${model} for ${dataset} with BO: ${bo}"
    if [ $running -ge $N_GPUS ]; then
      wait
      running=0
    fi
  done

done