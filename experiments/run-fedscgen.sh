#!/bin/bash
AVAILABLE_GPUS="${1:-0,1,2,3}"

declare -a TASK_QUEUE


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
    dataset="${DATASETS[$index]}"
    if [ "$inclusion" != "all" ] && { [ "$dataset" == "HumanDendriticCells" ] || [ "$dataset" == "CellLine" ]; }; then
        continue
    fi
    if [ "$dataset" == "HumanPancreas" ]; then
      batches="0,1,2,3,4"
      if [ "$inclusion" == "all" ]; then
        n_clients="5 4"
        batch_out="0 1"
      else
        n_clients="5"
        batch_out="0"
      fi
    elif [ "$dataset" == "CellLine" ]; then
      n_clients="3"
      batches="0,1,2"
      batch_out="0"
    else
      n_clients="2"
      batches="0,1"
      batch_out="0"
    fi
    combined=$([ "$inclusion" == "combined" ] && echo true || echo false)
    dropped=$([ "$inclusion" == "dropped" ] && echo true || echo false)
    task_name="${dataset}-${inclusion}"
    task="$task_name|${dataset}.h5ad|${DROPPED_CELLTYPES[$index]:-''}|$combined|$dropped|$batch_out|$n_clients|$batches|_GPU_"
    TASK_QUEUE+=("$task")
    break
  done
  break
done

chmod +x gpumaestro.sh
chmod +x fedscgen.sh
./gpumaestro.sh "$AVAILABLE_GPUS" "./fedscgen.sh" "${TASK_QUEUE[@]}"