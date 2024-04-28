#!/bin/bash

NUM_GPUS=3
GPU=0

# Making the scgen.sh executable
chmod +x scgen.sh


DATASETS=(HumanDendriticCells MouseCellAtlas HumanPancreas PBMC CellLine MouseRetina MouseBrain MouseHematopoieticStemProgenitorCells)
DROPPED_CELLTYPES=( ""
 "Epithelial,Dendritic,Smooth-muscle,NK"
  "stellate,endothelial,mesenchymal,macrophage,mast,epsilon,schwann,t_cell,MHC class II"
   "Plasmacytoid dendritic cell,Megakaryocyte,Hematopoietic stem cell"
   ""
   "ganglion,vascular_endothelium,horizontal,fibroblasts,microglia,pericytes,astrocytes"
   "Olfactory ensheathing cells,Choroid_plexus,Mitotic"
   "MPP,LTHSC,LMPP,Unsorted")

for index in "${!DATASETS[@]}"
do
  for inclusion in all dropped combined
  do
    dataset=${DATASETS[$index]}
    dropped_celltypes=${DROPPED_CELLTYPES[$index]}
    if [ $inclusion != "all" ]; then
      if [ "$dataset" == "HumanDendriticCells" ] || [ "$dataset" == "CellLine" ]; then
        continue
      fi
    fi
    combined=false
    dropped=false
    if [ $inclusion == "combined" ]; then
      combined=true
    elif [ $inclusion == "dropped" ]; then
      dropped=true
    else
      dropped_celltypes=""
    fi
    echo -e "\e[31mRunning scgen for $dataset with $inclusion with combined=$combined and dropped=$dropped and dropped_celltypes=$dropped_celltypes\e[0m"
    ./scgen.sh "$dataset.h5ad" "${dropped_celltypes}" $combined $dropped $GPU &
    GPU=$((GPU+1))
    if [ $GPU -eq $NUM_GPUS ]; then
      wait
      GPU=0
    fi

  done
done
wait