# config.sh

#SEEDS=(42 123 456 789 101112 2024 31415 2718 1618 9999)
#
#DATASETS=(HumanDendriticCells MouseCellAtlas HumanPancreas PBMC CellLine MouseRetina MouseBrain MouseHematopoieticStemProgenitorCells)
CELLTYPES=( ""
 "Epithelial,Dendritic,Smooth-muscle,NK"
  "stellate,endothelial,mesenchymal,macrophage,mast,epsilon,schwann,t_cell,MHC class II"
   "Plasmacytoid dendritic cell,Megakaryocyte,Hematopoietic stem cell"
   ""
   "ganglion,vascular_endothelium,horizontal,fibroblasts,microglia,pericytes,astrocytes"
   "Olfactory ensheathing cells,Choroid_plexus,Mitotic"
   "MPP,LTHSC,LMPP,Unsorted")

SEEDS=(123)

DATASETS=(MouseCellAtlas)

DROPPED_DATASETS=(MouseCellAtlas HumanPancreas PBMC MouseRetina MouseBrain MouseHematopoieticStemProgenitorCells)

DROPPED_CELLTYPES=(
  "Epithelial,Dendritic,Smooth-muscle,NK"
  "stellate,endothelial,mesenchymal,macrophage,mast,epsilon,schwann,t_cell,MHC class II"
  "Plasmacytoid dendritic cell,Megakaryocyte,Hematopoietic stem cell"
  "ganglion,vascular_endothelium,horizontal,fibroblasts,microglia,pericytes,astrocytes"
  "Olfactory ensheathing cells,Choroid_plexus,Mitotic"
  "MPP,LTHSC,LMPP,Unsorted"
)
