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
    combined=$([ "$inclusion" == "combined" ] && echo true || echo false)
    dropped=$([ "$inclusion" == "dropped" ] && echo true || echo false)
    n_clients=$([ "${DATASETS[$index]}" == "HumanPancreas" ] && echo "5 4" || echo "2")
    batches=$([ "${DATASETS[$index]}" == "HumanPancreas" ] && echo "0,1,2,3,4" || echo "0,1")
    batch_out=$([ "${DATASETS[$index]}" == "HumanPancreas" ] && echo "0 1" || echo "0")

    echo -e "\e[31mRunning fedscgen for ${DATASETS[$index]} with $inclusion with combined=$combined and dropped=$dropped and dropped_celltypes=${DROPPED_CELLTYPES[$index]} and n_clients=$n_clients and batches=$batches on GPU $GPU\e[0m"
    ./fedscgen.sh "${DATASETS[$index]}.h5ad" "${DROPPED_CELLTYPES[$index]}" $combined $dropped "$batch_out" "$n_clients" "$batches" "$GPU"&
    GPU=$((GPU+1))
    if [ $GPU -eq $NUM_GPUS ]; then
      wait
      GPU=0
    fi
  done
done
wait