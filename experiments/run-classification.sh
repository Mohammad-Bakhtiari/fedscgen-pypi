#!/bin/bash

# Making the classification.sh executable
chmod +x classification.sh

# 1: HumanDendriticCells
# Scenario: Only two classes, so no cell types to remove
## Model: mlp-norm
./classification.sh "all" "HumanDendriticCells" "1" "" "mlp-norm"
## Model: knn
./classification.sh "all" "HumanDendriticCells" "1" "" "knn"

## 2: MouseCellAtlas
## Model: mlp-norm
./classification.sh "all" "MouseCellAtlas" "1" "" "mlp-norm"

# Model: knn
./classification.sh "all" "MouseCellAtlas" "1" "" "knn"


## 4: HumanPancreas
## Model: mlp-norm
./classification.sh "all" "HumanPancreas" "4" "" "mlp-norm"

# Model: knn
./classification.sh "all" "HumanPancreas" "4" "" "knn"

# 5: PBMC
# Model: mlp-norm
./classification.sh "all" "PBMC" "1" "" "mlp-norm"

# Model: knn
./classification.sh "all" "PBMC" "1" "" "knn"


# 6: CellLine
# Model: mlp-norm
./classification.sh "all" "CellLine" "2" "" "mlp-norm"
# Model: knn
./classification.sh "all" "CellLine" "2" "" "knn"

# 7: MouseRetina
# Model: mlp-norm
./classification.sh "all" "MouseRetina" "1" "" "mlp-norm"

# Model: knn
./classification.sh "all" "MouseRetina" "1" "" "knn"

# 8: MouseBrain
# Model: mlp-norm
./classification.sh "all" "MouseBrain" "1" "" "mlp-norm"
# Model: knn
./classification.sh "all" "MouseBrain" "1" "" "knn"

## 10: MouseHematopoieticStemProgenitorCells
## Model: mlp-norm
./classification.sh "all" "MouseHematopoieticStemProgenitorCells" "1" "" "mlp-norm"
## Model: knn
./classification.sh "all" "MouseHematopoieticStemProgenitorCells" "1" "" "knn"
