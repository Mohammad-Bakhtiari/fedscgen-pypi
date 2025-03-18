#!/bin/bash
AVAILABLE_GPUS="${1:-0,1,2,3}"

declare -a TASK_QUEUE


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
  for inclusion in all dropped combined
  do
    dataset="${DATASETS[$index]}"
    combined=$([ "$inclusion" == "combined" ] && echo true || echo false)
    dropped=$([ "$inclusion" == "dropped" ] && echo true || echo false)
    n_clients=$([ "${dataset}" == "HumanPancreas" ] && echo "5 4" || echo "2")
    batches=$([ "${dataset}" == "HumanPancreas" ] && echo "0,1,2,3,4" || echo "0,1")
    batch_out=$([ "${dataset}" == "HumanPancreas" ] && echo "0 1" || echo "0")
    if [ "${dataset}" == "CellLine" ]; then
      n_clients="3"
      batches="0,1,2"
    fi
    task_name="${DATASETS[$index]}-${inclusion}"
    task="$task_name|${DATASETS[$index]}.h5ad|${DROPPED_CELLTYPES[$index]:-''}|$combined|$dropped|$batch_out|$n_clients|$batches|_GPU_"
    TASK_QUEUE+=("$task")
  done
done

chmod +x gpumaestro.sh
chmod +x fedscgen.sh
./gpumaestro.sh "$AVAILABLE_GPUS" "./fedscgen.sh" "${TASK_QUEUE[@]}"