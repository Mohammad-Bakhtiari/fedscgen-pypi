#!/bin/bash
NUM_GPUS="${2:-3}"

chmod +x fedscgen.sh

GPU=0
DATASETS=(MouseCellAtlas HumanPancreas PBMC MouseRetina MouseBrain MouseHematopoieticStemProgenitorCells)

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
  for inclusion in dropped combined
  do
    dataset=${DATASETS[$index]}
    dropped_celltypes=${DROPPED_CELLTYPES[$index]}
    combined=$([ "$inclusion" == "combined" ] && echo true || echo false)
    dropped=$([ "$inclusion" == "dropped" ] && echo true || echo false)
    n_clients="2"
    batches="0,1"
    batch_out="0"
    if [ "$dataset" == "HumanPancreas" ]; then
      n_clients="5 4"
      batches="0,1,2,3,4"
      batch_out="0 1"
    fi
    echo -e "\e[31mRunning fedscgen for $dataset with $inclusion with combined=$combined and dropped=$dropped and dropped_celltypes=$dropped_celltypes and n_clients=$n_clients and batches=$batches on GPU $GPU\e[0m"
    ./fedscgen.sh "$dataset.h5ad" "${dropped_celltypes}" $combined $dropped "$batch_out" "$n_clients" "$batches" "$GPU"&
    GPU=$((GPU+1))
    if [ $GPU -eq $NUM_GPUS ]; then
      wait
      GPU=0
    fi
  done
done
wait