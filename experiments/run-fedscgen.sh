#!/bin/bash

chmod +x fedscgen.sh
# HumanDendriticCells.h5ad
# It has only two classes, so Scenario 1 (dropping) and Scenario 2 (combining) does not work here!
# Scenario 3: Neither combining nor dropping any cell types
./fedscgen.sh "HumanDendriticCells.h5ad" "" false false "0" "2" "0,1" 0

# MouseCellAtlas.h5ad

# Scenario 1: Dropping cell types but not combining
# The cell types "Epithelial,Dendritic,Smooth-muscle,NK" will be dropped.
# Combining is set to false.
./fedscgen.sh "MouseCellAtlas.h5ad" "Epithelial,Dendritic,Smooth-muscle,NK" false true "0" "2" "0,1" 1

# Scenario 2: Combining cell types but not dropping
# The cell types "Epithelial,Dendritic,Smooth-muscle,NK" will be combined.
# Dropping is set to false.
./fedscgen.sh "MouseCellAtlas.h5ad" "Epithelial,Dendritic,Smooth-muscle,NK" true false "0" "2" "0,1" 1

# Scenario 3: Neither combining nor dropping any cell types
./fedscgen.sh "MouseCellAtlas.h5ad" "" false false "0" "2" "0,1" 0


# HumanPancreas.h5ad
# Scenario 1: Dropping cell types but not combining
./fedscgen.sh "HumanPancreas.h5ad" "stellate,endothelial,mesenchymal,macrophage,mast,epsilon,schwann,t_cell,MHC class II" false true "0 1 2 3" "5 4 3 2" "0,1,2,3,4" 1

# Scenario 2: Combining cell types but not dropping
./fedscgen.sh "HumanPancreas.h5ad" "stellate,endothelial,mesenchymal,macrophage,mast,epsilon,schwann,t_cell,MHC class II" true false "0 1 2 3" "5 4 3 2" "0,1,2,3,4" 1

# Scenario 3: Neither combining nor dropping any cell types
./fedscgen.sh "HumanPancreas.h5ad" "" false false "0 1 2 3" "5 4 3 2" "0,1,2,3,4" 1

# PBMC.h5ad
# Scenario 1: Dropping cell types but not combining
./fedscgen.sh "PBMC.h5ad" "Plasmacytoid dendritic cell,Megakaryocyte,Hematopoietic stem cell" false true "0" "2" "0,1" 1

# Scenario 2: Combining cell types but not dropping
./fedscgen.sh "PBMC.h5ad" "Plasmacytoid dendritic cell,Megakaryocyte,Hematopoietic stem cell" true false "0" "2" "0,1" 1

# Scenario 3: Neither combining nor dropping any cell types
./fedscgen.sh "PBMC.h5ad" "" false false "0" "2" "0,1" 1

# CellLine.h5ad
# This dataset doesn't have any specific cell types to combine or drop, so only Scenario 3 applies.
# Scenario 3: Neither combining nor dropping any cell types
./fedscgen.sh "CellLine.h5ad" "" false false "0 1" "3 2" "0,1,2" 1

# MouseRetina.h5ad
# Scenario 1: Dropping cell types but not combining
./fedscgen.sh "MouseRetina.h5ad" "ganglion,vascular_endothelium,horizontal,fibroblasts,microglia,pericytes,astrocytes" false true "0" "2" "0,1" 1

# Scenario 2: Combining cell types but not dropping
./fedscgen.sh "MouseRetina.h5ad" "ganglion,vascular_endothelium,horizontal,fibroblasts,microglia,pericytes,astrocytes" true false "0" "2" "0,1" 1

# Scenario 3: Neither combining nor dropping any cell types
./fedscgen.sh "MouseRetina.h5ad" "" false false "0" "2" "0,1" 1

# MouseBrain.h5ad
# Scenario 1: Dropping cell types but not combining
./fedscgen.sh "MouseBrain.h5ad" "Olfactory ensheathing cells,Choroid_plexus,Mitotic" false true "0" "2" "0,1" 1

# Scenario 2: Combining cell types but not dropping
./fedscgen.sh "MouseBrain.h5ad" "Olfactory ensheathing cells,Choroid_plexus,Mitotic" true false "0" "2" "0,1" 1

# Scenario 3: Neither combining nor dropping any cell types
./fedscgen.sh "MouseBrain.h5ad" "" false false "0" "2" "0,1" 1

# MouseHematopoieticStemProgenitorCells.h5ad
# Scenario 1: Dropping cell types but not combining
./fedscgen.sh "MouseHematopoieticStemProgenitorCells.h5ad" "MPP,LTHSC,LMPP,Unsorted" false true "0" "2" "0,1" 1

# Scenario 2: Combining cell types but not dropping
./fedscgen.sh "MouseHematopoieticStemProgenitorCells.h5ad" "MPP,LTHSC,LMPP,Unsorted" true false "0" "2" "0.1" 1

# Scenario 3: Neither combining nor dropping any cell types
./fedscgen.sh "MouseHematopoieticStemProgenitorCells.h5ad" "" false false "0" "2" "0,1" 1