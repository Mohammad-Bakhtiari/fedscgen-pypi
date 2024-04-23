import __init__
import os
import random

random.seed(42)

from fedscgen.utils import to_csv, plot_metrics, encode_labels, normalize_data, classify_celltypes, custom_kfold,\
    remove_cell_types, combine_cell_types, get_cuda_device
from sklearn.model_selection import StratifiedKFold
import pandas as pd
import numpy as np
import argparse
import anndata
import torch
torch.multiprocessing.set_sharing_strategy('file_system')


def cross_validate(x, y, epochs, batch_size, lr, output_dir, model, hidden_size, batches, batch_out, init_model=None, device='cpu'):
    y = encode_labels(y)
    n_classes = len(np.unique(y))

    if batch_out > 0:
        kf = custom_kfold(batches, batch_out)
    else:
        kf = StratifiedKFold(n_splits=5).split(x, y)

    metrics_dfs = []

    for fold, (train_ind, test_ind) in enumerate(kf, 1):
        x_train, y_train, x_test, y_test = x[train_ind], y[train_ind], x[test_ind], y[test_ind]
        metrics = classify_celltypes(x_train, y_train, x_test, y_test, epochs, lr, batch_size, n_classes, init_model,
                                     model, hidden_size, device)

        filename = os.path.join(output_dir, f"classification_acc_{fold}.csv")
        df = to_csv(*metrics, filename)
        df['Fold'] = fold  # add 'Fold' column to keep track of the fold
        metrics_dfs.append(df)

    # concatenate all dataframes and compute the mean metrics
    metrics_all = pd.concat(metrics_dfs)
    metrics_summary = metrics_all.groupby(['Epoch', 'Fold']).mean().reset_index()
    plot_metrics(metrics_summary[["Epoch", "AUC", "Accuracy"]], plt_name=os.path.join(output_dir, f"mean_metrics.png"))
    if model.lower() == "knn":
        plot_metrics(metrics_summary[['Epoch', 'Train Loss', 'Test Loss']],
                     plt_name=os.path.join(output_dir, f"mean_loss.png"))


if __name__ == '__main__':
    """
        This script should be called for only one H5AD file. 
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("--init_model_path", type=str, default="./model_repo/mlp.pth")
    parser.add_argument("--adata", type=str, default="./datasets/reproduce/pancreas-cleaned.h5ad")
    parser.add_argument("--output", type=str, default="./results/classification/centralized/pancreas/uncorrected")
    parser.add_argument("--epoch", type=int, default=50)
    parser.add_argument("--cell_key", type=str, default="celltype")
    parser.add_argument("--batch_key", type=str, default="sample")
    parser.add_argument("--lr", type=float, default=0.001)
    parser.add_argument("--batch_size", type=int, default=32)
    # Define cell types to be removed
    # remove_cell_types = "not applicable,unclassified endocrine,unclassified,unclear"
    parser.add_argument("--remove_cell_types", type=str, default="")
    parser.add_argument("--model", type=str, default="mlp-norm", choices=["mlp", "mlp-norm", "knn", "kmeans"])
    parser.add_argument("--hidden_size", type=str, default="800,800")
    parser.add_argument("--overwrite", action='store_true', default=False)
    parser.add_argument('--norm_method', type=str, default="z_score",
                        help='Normalization method: "log", "quantile", "z_score", "min_max", or "tpm".')
    parser.add_argument("--batch_out", type=int, default=0)
    parser.add_argument("--combine", action='store_true', default=False)
    parser.add_argument("--gpu", type=int, default=0)

    args = parser.parse_args()
    args.device = get_cuda_device(args.gpu)
    if args.model == "mlp-norm":
        args.init_model_path = args.init_model_path.replace("mlp", f"mlp-{args.hidden_size}")
    args.hidden_size = [int(num) for num in args.hidden_size.split(",")]
    args.remove_cell_types = [c.strip() for c in args.remove_cell_types.strip().split(",")]

    os.makedirs(args.output, exist_ok=True)
    adata = anndata.read_h5ad(args.adata)
    if args.combine:
        adata = combine_cell_types(adata, args.remove_cell_types, args.cell_key)
    else:
        adata = remove_cell_types(adata, args.remove_cell_types, args.cell_key)

    if args.model != "mlp":
        args.init_model_path = args.init_model_path.replace("mlp", args.model)

    x = adata.obsm["pca_20"]
    if len(args.norm_method) > 0:
        x = normalize_data(x, args.norm_method)
    y = adata.obs[args.cell_key]
    cross_validate(x=x,
                   y=y,
                   epochs=args.epoch,
                   batch_size=args.batch_size,
                   lr=args.lr,
                   output_dir=args.output,
                   init_model=args.init_model_path,
                   model=args.model,
                   hidden_size=args.hidden_size,
                   batches=adata.obs[args.batch_key],
                   batch_out=args.batch_out,
                   device=args.device
                   )
