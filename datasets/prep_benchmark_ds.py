import argparse
from utils import prep

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--root_path", type=str, required=True)
    parser.add_argument("--dataset", type=str, required=True, choices=[f"dataset{i}" for i in [1, 2, 4, 5, 6, 7, 8, 9, 10]])
    args = parser.parse_args()
    datasets = {
        "dataset1": {"name": "HumanDendriticCells",
                     "adata_filename": f"HumanDendriticCells.h5ad",
                     "expression": "dataset1_sm_uc3.txt",
                     "annotation": "sample_sm_uc3.txt",
                     "batch_key": "batch",
                     "cell_key": "celltype",
                     },
        "dataset2": {"name": "MouseCellAtlas",
                     "adata_filename": f"MouseCellAtlas.h5ad",
                     "expression": "filtered_total_batch1_seqwell_batch2_10x.txt",
                     "annotation": "filtered_total_sample_ext_organ_celltype_batch.txt",
                     "batch_key": "batch",
                     "cell_key": "ct"
                     },
        # We skip the third dataset
        # "dataset3": {"name": "simulation",
        #              "adata_filename": f"simulation.h5ad",
        #              "expression": "dataset1_sm_uc3.txt",
        #              "annotation": "sample_sm_uc3.txt",
        #              "batch_key": "batch",
        #              "cell_key": "celltype"
        #              },
        "dataset4": {"name": "HumanPancreas",
                     "adata_filename": f"HumanPancreas.h5ad",
                     "expression": "myData_pancreatic_5batches.txt",
                     "annotation": "mySample_pancreatic_5batches.txt",
                     "batch_key": "batchlb",
                     "cell_key": "celltype",
                     },
        "dataset5": {"name": "PBMC",
                     "adata_filename": f"PBMC.h5ad",
                     "expression": "b1_exprs.txt,b2_exprs.txt",
                     "annotation": "b1_celltype.txt,b2_celltype.txt",
                     "batch_key": "batch",
                     "cell_key": "CellType",
                     },
        "dataset6": {"name": "CellLine",
                     "adata_filename": f"CellLine.h5ad",
                     "expression": "b1_exprs.txt,b2_exprs.txt,b3_exprs.txt",
                     "annotation": "b1_celltype.txt,b2_celltype.txt,b3_celltype.txt",
                     "batch_key": "None",
                     "cell_key": "CellType",
                     },
        "dataset7": {"name": "MouseRetina",
                     "adata_filename": f"MouseRetina.h5ad",
                     "expression": "b1_exprs.txt,b2_exprs.txt",
                     "annotation": "b1_celltype.txt,b2_celltype.txt",
                     "batch_key": "None",
                     "cell_key": "CellType",
                     },
        "dataset8": {"name": "MouseBrain",
                     "adata_filename": f"MouseBrain.h5ad",
                     "expression": "dropviz_and_nuclei_combined_filtered_UMI.mtx",
                     "annotation": "dropviz_and_nuclei_combined_filtered_cell_info.txt",
                     "batch_key": "batch",
                     "cell_key": "cell_type",
                     },
        # We skip the ninth dataset
        # "dataset9": {"name": "HumanCellAtlas",
        #              "adata_filename": f"HumanCellAtlas.h5ad",
        #              "expression": "HCA_genes_cells_filtered_filtered_UMI.mtx",
        #              "annotation": "HCA_genes_cells_filtered_filtered_cell_info.txt",
        #              "batch_key": "batch",
        #              "cell_key": "cell_type",
        #              },
        "dataset10": {"name": "MouseHematopoieticStemProgenitorCells",
                      "adata_filename": f"MouseHematopoieticStemProgenitorCells.h5ad",
                      "expression": "b1_exprs.txt,b2_exprs.txt",
                      "annotation": "b1_celltype.txt,b2_celltype.txt",
                      "batch_key": "None",
                      "cell_key": "CellType",
                      }
    }
    # dataset_preprocessors = getattr(benchmark_dataset_preprocessors, f"{args.dataset}_prep")
    prep(path=args.root_path, **datasets[args.dataset])