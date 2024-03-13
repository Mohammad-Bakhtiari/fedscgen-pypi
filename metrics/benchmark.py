import numpy as np
import pandas as pd
import anndata
import argparse
import os
from scarches.metrics import nmi, entropy_batch_mixing, asw
from utils import (graph_connectivity_score, isolated_label_f1_score, ari_score, bar_plot, plot_metrics_with_circles,
                   compute_ils, knn_accuracy, DATASETS, bar_plot_subplot)


def calculate_and_plot_metrics(adata_dict, batch_key, cell_key, plot_name, overwrite=False, n_components=50):
    metric_file = plot_name.replace(".png", ".csv")
    if os.path.exists(metric_file) and not overwrite:
        df = pd.read_csv(metric_file, sep=',')
    else:
        all_metrics = []
        for key, adata in adata_dict.items():
            latent_adata = anndata.AnnData(adata.obsm[f'pca_{n_components}'])
            latent_adata.obs = adata.obs
            all_metrics.append(benchmark(adata, latent_adata, batch_key, cell_key, key))
    plot_and_save(all_metrics, plot_name)
    #     df = pd.DataFrame(all_metrics)
    #     df.to_csv(metric_file, sep=",", index=False)
    #     print(df)
    # bar_plot(df, plot_name.replace(".", "-bar."))
    # plot_metrics_with_circles(df, plot_name.replace(".", "-circle."))


def benchmark(adata, latent_adata, batch_key, cell_key, dataset_name):
    nmi_value = nmi(latent_adata, cell_key)
    gc_value = graph_connectivity_score(latent_adata)
    f1_value = isolated_label_f1_score(latent_adata, cell_key, batch_key)
    ari_value = ari_score(adata, cell_key)
    asw_b, asw_c = asw(latent_adata, cell_key, batch_key)
    ebm_value = entropy_batch_mixing(adata, batch_key, n_neighbors=15)
    ils_value = compute_ils(adata, cell_key, batch_key)
    knn_acc_value = knn_accuracy(latent_adata, cell_key)
    common_metrics = {
        'Dataset': dataset_name,
        'NMI': nmi_value,
        'GC': gc_value,
        'ILF1': f1_value,
        'ARI': ari_value,
        'EBM': ebm_value,
        'ILS': ils_value,
        'KNN Acc': knn_acc_value,
        'ASW_B': asw_b,
        'ASW_C': asw_c
    }
    return common_metrics


def benchmark_nmi_ari_ebm(adata, latent_adata, batch_key, cell_key, epoch, round, approach):
    nmi_value = nmi(latent_adata, cell_key)
    ari_value = ari_score(adata, cell_key)
    ebm_value = entropy_batch_mixing(adata, batch_key, n_neighbors=15)
    common_metrics = {
        'Approach': approach,
        'Epoch': epoch,
        'Round': round,
        'NMI': nmi_value,
        'ARI': ari_value,
        'EBM': ebm_value,
    }
    return common_metrics


def benchmark_snapshots(data_dir, n_rounds, n_components, batch_key, cell_key):
    """

    Parameters
    ----------
    data_dir: str
        Full path to the main data directory
    n_rounds:
        number of rounds

    Returns
    -------

    """
    all_metrics = []
    for r in range(1, n_rounds + 1):
        adata = anndata.read_h5ad(os.path.join(data_dir, f"FedScGen-C5-{r}.h5ad"))
        latent_adata = anndata.AnnData(adata.obsm[f'pca_{n_components}'])
        latent_adata.obs = adata.obs
        all_metrics.append(benchmark(adata, latent_adata, batch_key, cell_key, str(r)))
    adata = anndata.read_h5ad(os.path.join(data_dir, f"scgen.h5ad"))
    latent_adata = anndata.AnnData(adata.obsm[f'pca_{n_components}'])
    latent_adata.obs = adata.obs
    all_metrics.append(benchmark(adata, latent_adata, batch_key, cell_key, "ScGen"))
    plot_and_save(all_metrics, os.path.join(data_dir, "metrics.png"))


def plot_and_save(all_metrics, plot_name):
    df = pd.DataFrame(all_metrics)
    df.to_csv(plot_name.replace(".png", ".csv"), sep=",", index=False)
    bar_plot(df, plot_name.replace(".", "-bar."))
    plot_metrics_with_circles(df, plot_name.replace(".", "-circle."))


def benchmark_all_datasets(data_dir: str, datasets: list, variant: str, n_components: int):
    for ds in datasets:
        dataset_path = os.path.join(data_dir, ds, variant)
        if not os.path.exists(dataset_path):
            raise FileNotFoundError(f"{dataset_path} Dose not exist")

        adata_dict = {
            'Raw': anndata.read_h5ad(os.path.join(dataset_path, 'Raw.h5ad')),
            'ScGen': anndata.read_h5ad(os.path.join(dataset_path, 'ScGen.h5ad')),
            'FedScGen': anndata.read_h5ad(os.path.join(dataset_path, 'FedScGen.h5ad'))
        }
        plot_name_variant = os.path.join(dataset_path, "metrics.png")
        calculate_and_plot_metrics(adata_dict, "batch", "cell_type", plot_name_variant, overwrite=True,
                                   n_components=n_components)


def benchmark_batch_out(data_dir, n_batches, n_components, batch_key, cell_key):
    for b in range(n_batches):
        all_metrics = []
        for model in ["ScGen", "FedScGen"]:
            dataset_path = os.path.join(data_dir, str(b))
            adata = anndata.read_h5ad(os.path.join(dataset_path, f"{model}.h5ad"))
            latent_adata = anndata.AnnData(adata.obsm[f'pca_{n_components}'])
            latent_adata.obs = adata.obs
            all_metrics.append(benchmark(adata, latent_adata, batch_key, cell_key, model))
        plot_name = os.path.join(dataset_path, "metrics.png")
        plot_and_save(all_metrics, plot_name)


def benchmark_tuning(data_dir, n_components, batch_key, cell_key):
    """
    benchmark all files in the tuning directory
    Returns
    -------

    """
    h5ad_files = {}
    for epoch in range(1, 11):
        files = {}
        for round in range(1, 11):
            adata = anndata.read_h5ad(os.path.join(data_dir, f"E{epoch}", f"corrected_{round}.h5ad"))
            files[round] = adata
        files = sorted(files.items(), key=lambda x: x[0])
        h5ad_files[epoch] = files
    all_metrics = []
    sorted_files = sorted(h5ad_files.items(), key=lambda x: x[0])
    for epoch, file in sorted_files:
        for round, adata in file:
            latent_adata = anndata.AnnData(adata.obsm[f'pca_{n_components}'])
            latent_adata.obs = adata.obs
            all_metrics.append(
                benchmark_nmi_ari_ebm(adata, latent_adata, batch_key, cell_key, epoch, round, 'FedscGen'))
    # add scGen
    adata = anndata.read_h5ad(os.path.join(data_dir, "scGen.h5ad"))
    latent_adata = anndata.AnnData(adata.obsm[f'pca_{n_components}'])
    latent_adata.obs = adata.obs
    all_metrics.append(benchmark_nmi_ari_ebm(adata, latent_adata, batch_key, cell_key, 0, 0, f"scGen"))
    df = pd.DataFrame(all_metrics)
    plot_name = os.path.join(data_dir, "metrics.png")
    df.to_csv(os.path.join(data_dir, "metrics.csv"), sep=",", index=False)
    bar_plot_subplot(df, plot_name)


def load_metrics_and_plot(df_path, plot_name):
    df = pd.read_csv(df_path)
    bar_plot(df, plot_name.replace(".", "-bar."))
    df.drop(columns=["ILS", "ILF1"], inplace=True)
    plot_metrics_with_circles(df, plot_name.replace(".", "-circle."))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate and Plot NMI for a set of adata files.')

    # Add arguments
    parser.add_argument("--data_dir", type=str, help="path to the main data directory",
                        default="/home/bba1658/FedSC/results/scgen/batchout")
    parser.add_argument('--n_components', type=int, default=20, help="Number of components for PCA")
    parser.add_argument("--inclusion", type=str, default="all", choices=["all", "combined", "dropped"])
    parser.add_argument('--cell_key', help='Cell key name.', default="cell_type")
    parser.add_argument('--batch_key', help='Batch key name.', default="batch")
    parser.add_argument("--scenarios", type=str, default="all",
                        choices=["datasets", "batch-out", "snapshots", "tuning"])
    parser.add_argument("--n_rounds", type=int, default=10, help="Number of rounds for snapshots")
    parser.add_argument("--n_batches", type=int, default=10, help="Number of batches for batch-out")
    parser.add_argument("--plot_only", action="store_true", default=False, help="Plot the metrics")
    args = parser.parse_args()

    if args.plot_only:
        if args.scenarios == "tuning":
            df = pd.read_csv(os.path.join(args.data_dir, "metrics.csv"))
            bar_plot_subplot(df, os.path.join(args.data_dir, "metrics.png"))
        else:
            load_metrics_and_plot(os.path.join(args.data_dir, "metrics.csv"),
                                  os.path.join(args.data_dir, "metrics.png"))
    else:
        if args.scenarios == "datasets":
            benchmark_all_datasets(args.data_dir, DATASETS, args.inclusion, args.n_components)
        elif args.scenarios == "batch-out":
            benchmark_batch_out(args.data_dir, args.n_batches, args.n_components, args.batch_key, args.cell_key)
        elif args.scenarios == "tuning":
            benchmark_tuning(args.data_dir, args.n_components, args.batch_key, args.cell_key)
        else:
            benchmark_snapshots(args.data_dir, args.n_rounds, args.n_components, args.batch_key, args.cell_key)
