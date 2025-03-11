#!/bin/bash

NUM_GPUS="${1:-3}"
GPU=0

# Making the scgen.sh executable
chmod +x scgen.sh


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
    combined=$([ "$inclusion" == "combined" ] && echo true || echo false)
    dropped=$([ "$inclusion" == "dropped" ] && echo true || echo false)
    echo -e "\e[31mRunning scgen for ${DATASETS[$index]} with $inclusion with combined=$combined and dropped=$dropped and dropped_celltypes=${DROPPED_CELLTYPES[$index]}\e[0m"
    ./scgen.sh "${DATASETS[$index]}.h5ad" "${DROPPED_CELLTYPES[$index]}" "$combined" "$dropped" $GPU &
    GPU=$((GPU+1))
    if [ $GPU -eq $NUM_GPUS ]; then
      wait
      GPU=0
    fi

  done
done
wait