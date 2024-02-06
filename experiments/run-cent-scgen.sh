#!/bin/bash

# call_general_script.sh
# This script acts as a wrapper to run_script.sh, allowing for multiple runs with different parameters.

# Making the cent-scgen.sh executable
chmod +x cent-scgen.sh

# Dataset1: HumanDendriticCells.h5ad
# It has only two classes, so Scenario 1 (dropping) and Scenario 2 (combining) does not work here!
# Scenario 3: Neither combining nor dropping any cell types for dataset1
./cent-scgen.sh "dataset1" "HumanDendriticCells.h5ad" "" false false

# Dataset2: MouseCellAtlas.h5ad

# Scenario 1: Dropping cell types but not combining for dataset2
# The cell types "Epithelial,Dendritic,Smooth-muscle,NK" will be dropped.
# Combining is set to false.
./cent-scgen.sh "dataset2" "MouseCellAtlas.h5ad" "Epithelial,Dendritic,Smooth-muscle,NK" false true

# Scenario 2: Combining cell types but not dropping for dataset2
# The cell types "Epithelial,Dendritic,Smooth-muscle,NK" will be combined.
# Dropping is set to false.
./cent-scgen.sh "dataset2" "MouseCellAtlas.h5ad" "Epithelial,Dendritic,Smooth-muscle,NK" true false

# Scenario 3: Neither combining nor dropping any cell types for dataset2
./cent-scgen.sh "dataset2" "MouseCellAtlas.h5ad" "" false false


# Dataset4: HumanPancreas.h5ad
# Scenario 1: Dropping cell types but not combining for dataset4
./cent-scgen.sh "dataset4" "HumanPancreas.h5ad" "stellate,endothelial,mesenchymal,macrophage,mast,epsilon,schwann,t_cell,MHC class II" false true

# Scenario 2: Combining cell types but not dropping for dataset4
./cent-scgen.sh "dataset4" "HumanPancreas.h5ad" "stellate,endothelial,mesenchymal,macrophage,mast,epsilon,schwann,t_cell,MHC class II" true false

# Scenario 3: Neither combining nor dropping any cell types for dataset4
./cent-scgen.sh "dataset4" "HumanPancreas.h5ad" "" false false

# Dataset5: PBMC.h5ad
# Scenario 1: Dropping cell types but not combining for dataset5
./cent-scgen.sh "dataset5" "PBMC.h5ad" "Plasmacytoid dendritic cell,Megakaryocyte,Hematopoietic stem cell" false true

# Scenario 2: Combining cell types but not dropping for dataset5
./cent-scgen.sh "dataset5" "PBMC.h5ad" "Plasmacytoid dendritic cell,Megakaryocyte,Hematopoietic stem cell" true false

# Scenario 3: Neither combining nor dropping any cell types for dataset5
./cent-scgen.sh "dataset5" "PBMC.h5ad" "" false false

# Dataset6: CellLine.h5ad
# This dataset doesn't have any specific cell types to combine or drop, so only Scenario 3 applies.
# Scenario 3: Neither combining nor dropping any cell types for dataset6
./cent-scgen.sh "dataset6" "CellLine.h5ad" "" false false

# Dataset7: MouseRetina.h5ad
# Scenario 1: Dropping cell types but not combining for dataset7
./cent-scgen.sh "dataset7" "MouseRetina.h5ad" "ganglion,vascular_endothelium,horizontal,fibroblasts,microglia,pericytes,astrocytes" false true

# Scenario 2: Combining cell types but not dropping for dataset7
./cent-scgen.sh "dataset7" "MouseRetina.h5ad" "ganglion,vascular_endothelium,horizontal,fibroblasts,microglia,pericytes,astrocytes" true false

# Scenario 3: Neither combining nor dropping any cell types for dataset7
./cent-scgen.sh "dataset7" "MouseRetina.h5ad" "" false false

# Dataset8: MouseBrain.h5ad
# Scenario 1: Dropping cell types but not combining for dataset8
./cent-scgen.sh "dataset8" "MouseBrain.h5ad" "Olfactory ensheathing cells,Choroid_plexus,Mitotic" false true

# Scenario 2: Combining cell types but not dropping for dataset8
./cent-scgen.sh "dataset8" "MouseBrain.h5ad" "Olfactory ensheathing cells,Choroid_plexus,Mitotic" true false

# Scenario 3: Neither combining nor dropping any cell types for dataset8
./cent-scgen.sh "dataset8" "MouseBrain.h5ad" "" false false

# Dataset10: MouseHematopoieticStemProgenitorCells.h5ad
# Scenario 1: Dropping cell types but not combining for dataset10
./cent-scgen.sh "dataset10" "MouseHematopoieticStemProgenitorCells.h5ad" "MPP,LTHSC,LMPP,Unsorted" false true

# Scenario 2: Combining cell types but not dropping for dataset10
./cent-scgen.sh "dataset10" "MouseHematopoieticStemProgenitorCells.h5ad" "MPP,LTHSC,LMPP,Unsorted" true false

# Scenario 3: Neither combining nor dropping any cell types for dataset10
./cent-scgen.sh "dataset10" "MouseHematopoieticStemProgenitorCells.h5ad" "" false false

