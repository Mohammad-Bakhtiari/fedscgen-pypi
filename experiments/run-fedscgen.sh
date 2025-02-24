#!/bin/bash
SMPC="${1:-false}"
NUM_GPUS="${2:-3}"

chmod +x fedscgen.sh

GPU=0
DATASETS=(HumanDendriticCells MouseCellAtlas HumanPancreas PBMC CellLine MouseRetina MouseBrain MouseHematopoieticStemProgenitorCells)
DATASETS=(CellLine)
DROPPED_CELLTYPES=( ""
 "Epithelial,Dendritic,Smooth-muscle,NK"
  "stellate,endothelial,mesenchymal,macrophage,mast,epsilon,schwann,t_cell,MHC class II"
   "Plasmacytoid dendritic cell,Megakaryocyte,Hematopoietic stem cell"
   ""
   "ganglion,vascular_endothelium,horizontal,fibroblasts,microglia,pericytes,astrocytes"
   "Olfactory ensheathing cells,Choroid_plexus,Mitotic"
   "MPP,LTHSC,LMPP,Unsorted")
DROPPED_CELLTYPES=("")

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
    n_clients="2"
    batches="0,1"
    batch_out="0"
    if [ "$dataset" == "HumanPancreas" ]; then
      n_clients="5 4"
      batches="0,1,2,3,4"
      batch_out="0 1"
    elif [ "$dataset" == "CellLine" ]; then
      n_clients="3"
      batches="0,1,2"

    fi
    echo -e "\e[31mRunning fedscgen for $dataset with $inclusion with combined=$combined and dropped=$dropped and dropped_celltypes=$dropped_celltypes and n_clients=$n_clients and batches=$batches on GPU $GPU\e[0m"
    ./fedscgen.sh "$dataset.h5ad" "${dropped_celltypes}" $combined $dropped "$batch_out" "$n_clients" "$batches" "$GPU" "$SMPC" &
    GPU=$((GPU+1))
    if [ $GPU -eq $NUM_GPUS ]; then
      wait
      GPU=0
    fi
  done
done
wait