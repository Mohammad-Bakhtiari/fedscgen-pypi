# config.sh

AVAILABLE_GPUS="${1:-0,1,2,3}"

DROPPED_DATASETS=(MouseCellAtlas HumanPancreas PBMC MouseRetina MouseBrain MouseHematopoieticStemProgenitorCells)

DROPPED_CELLTYPES=(
  "Epithelial,Dendritic,Smooth-muscle,NK"
  "stellate,endothelial,mesenchymal,macrophage,mast,epsilon,schwann,t_cell,MHC class II"
  "Plasmacytoid dendritic cell,Megakaryocyte,Hematopoietic stem cell"
  "ganglion,vascular_endothelium,horizontal,fibroblasts,microglia,pericytes,astrocytes"
  "Olfactory ensheathing cells,Choroid_plexus,Mitotic"
  "MPP,LTHSC,LMPP,Unsorted"
)
