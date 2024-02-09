#!/bin/bash

# call_general_script.sh
# This script acts as a wrapper to run_script.sh, allowing for multiple runs with different parameters.

# Making the scgen.sh executable
chmod +x scgen.sh

# HumanDendriticCells.h5ad
# It has only two classes, so Scenario 1 (dropping) and Scenario 2 (combining) does not work here!
# Scenario 3: Neither combining nor dropping any cell types
./scgen.sh "HumanDendriticCells.h5ad" "" false false

# MouseCellAtlas.h5ad

# Scenario 1: Dropping cell types but not combining
# The cell types "Epithelial,Dendritic,Smooth-muscle,NK" will be dropped.
# Combining is set to false.
./scgen.sh "MouseCellAtlas.h5ad" "Epithelial,Dendritic,Smooth-muscle,NK" false true

# Scenario 2: Combining cell types but not dropping
# The cell types "Epithelial,Dendritic,Smooth-muscle,NK" will be combined.
# Dropping is set to false.
./scgen.sh "MouseCellAtlas.h5ad" "Epithelial,Dendritic,Smooth-muscle,NK" true false

# Scenario 3: Neither combining nor dropping any cell types
./scgen.sh "MouseCellAtlas.h5ad" "" false false


# HumanPancreas.h5ad
# Scenario 1: Dropping cell types but not combining
./scgen.sh "HumanPancreas.h5ad" "stellate,endothelial,mesenchymal,macrophage,mast,epsilon,schwann,t_cell,MHC class II" false true

# Scenario 2: Combining cell types but not dropping
./scgen.sh "HumanPancreas.h5ad" "stellate,endothelial,mesenchymal,macrophage,mast,epsilon,schwann,t_cell,MHC class II" true false

# Scenario 3: Neither combining nor dropping any cell types
./scgen.sh "HumanPancreas.h5ad" "" false false

# PBMC.h5ad
# Scenario 1: Dropping cell types but not combining
./scgen.sh "PBMC.h5ad" "Plasmacytoid dendritic cell,Megakaryocyte,Hematopoietic stem cell" false true

# Scenario 2: Combining cell types but not dropping
./scgen.sh "PBMC.h5ad" "Plasmacytoid dendritic cell,Megakaryocyte,Hematopoietic stem cell" true false

# Scenario 3: Neither combining nor dropping any cell types
./scgen.sh "PBMC.h5ad" "" false false

# CellLine.h5ad
# This dataset doesn't have any specific cell types to combine or drop, so only Scenario 3 applies.
# Scenario 3: Neither combining nor dropping any cell types
./scgen.sh "CellLine.h5ad" "" false false

# MouseRetina.h5ad
# Scenario 1: Dropping cell types but not combining
./scgen.sh "MouseRetina.h5ad" "ganglion,vascular_endothelium,horizontal,fibroblasts,microglia,pericytes,astrocytes" false true

# Scenario 2: Combining cell types but not dropping
./scgen.sh "MouseRetina.h5ad" "ganglion,vascular_endothelium,horizontal,fibroblasts,microglia,pericytes,astrocytes" true false

# Scenario 3: Neither combining nor dropping any cell types
./scgen.sh "MouseRetina.h5ad" "" false false

# MouseBrain.h5ad
# Scenario 1: Dropping cell types but not combining
./scgen.sh "MouseBrain.h5ad" "Olfactory ensheathing cells,Choroid_plexus,Mitotic" false true

# Scenario 2: Combining cell types but not dropping
./scgen.sh "MouseBrain.h5ad" "Olfactory ensheathing cells,Choroid_plexus,Mitotic" true false

# Scenario 3: Neither combining nor dropping any cell types
./scgen.sh "MouseBrain.h5ad" "" false false

# MouseHematopoieticStemProgenitorCells.h5ad
# Scenario 1: Dropping cell types but not combining
./scgen.sh "MouseHematopoieticStemProgenitorCells.h5ad" "MPP,LTHSC,LMPP,Unsorted" false true

# Scenario 2: Combining cell types but not dropping
./scgen.sh "MouseHematopoieticStemProgenitorCells.h5ad" "MPP,LTHSC,LMPP,Unsorted" true false

# Scenario 3: Neither combining nor dropping any cell types
./scgen.sh "MouseHematopoieticStemProgenitorCells.h5ad" "" false false