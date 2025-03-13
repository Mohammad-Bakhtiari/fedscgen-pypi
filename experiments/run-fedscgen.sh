#!/bin/bash
AVAILABLE_GPUS="${1:-0,1,2,3}"
SMPC="${2:-"true"}"

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
  for inclusion in dropped combined
  do
    combined=$([ "$inclusion" == "combined" ] && echo true || echo false)
    dropped=$([ "$inclusion" == "dropped" ] && echo true || echo false)
    n_clients=$([ "${DATASETS[$index]}" == "HumanPancreas" ] && echo "5 4" || echo "2")
    batches=$([ "${DATASETS[$index]}" == "HumanPancreas" ] && echo "0,1,2,3,4" || echo "0,1")
    batch_out=$([ "${DATASETS[$index]}" == "HumanPancreas" ] && echo "0 1" || echo "0")
    task_name="${DATASETS[$index]}-${inclusion}"
    task="$task_name|${DATASETS[$index]}.h5ad|${DROPPED_CELLTYPES[$index]:-''}|$combined|$dropped|$batch_out|$n_clients|$batches|_GPU_"
    TASK_QUEUE+=("$task")
  done
done

script_name="fedscgen.sh"
[ "$SMPC" == "true" ] && script_name="fedscgen-smpc.sh"
chmod +x gpumaestro.sh
chmod +x $script_name
./gpumaestro.sh "$AVAILABLE_GPUS" "./${script_name}" "${TASK_QUEUE[@]}"