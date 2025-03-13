#!/bin/bash

AVAILABLE_GPUS="${1:-0,1,2,3}"

declare -a TASK_QUEUE

DATASETS=(MouseCellAtlas HumanPancreas PBMC MouseRetina MouseBrain MouseHematopoieticStemProgenitorCells)
DROPPED_CELLTYPES=(""
  "stellate,endothelial,mesenchymal,macrophage,mast,epsilon,schwann,t_cell,MHC class II"
   "Plasmacytoid dendritic cell,Megakaryocyte,Hematopoietic stem cell"
   "ganglion,vascular_endothelium,horizontal,fibroblasts,microglia,pericytes,astrocytes"
   "Olfactory ensheathing cells,Choroid_plexus,Mitotic"
   "MPP,LTHSC,LMPP,Unsorted")

for index in "${!DATASETS[@]}"
do
  for inclusion in dropped combined
  do
    combined=$([ "$inclusion" == "combined" ] && echo true || echo false)
    dropped=$([ "$inclusion" == "dropped" ] && echo true || echo false)
    task_name="${DATASETS[$index]}-${inclusion}"
    task="$task_name|${DATASETS[$index]}.h5ad|${DROPPED_CELLTYPES[$index]:-}|$combined|$dropped|_GPU_"
    TASK_QUEUE+=("$task")


  done
done
wait

chmod +x gpumaestro.sh
chmod +x scgen.sh
./gpumaestro.sh "$AVAILABLE_GPUS" "./scgen.sh" "${TASK_QUEUE[@]}"