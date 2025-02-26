import numpy as np
import pandas as pd
import anndata
import argparse
import os
from scarches.metrics import nmi, entropy_batch_mixing, asw
from scipy.stats import wilcoxon
from utils import (graph_connectivity_score, isolated_label_f1_score, ari_score, bar_plot, plot_metrics_with_circles,
                   compute_ils, knn_accuracy, DATASETS, bar_plot_subplot)
import sys
from pathlib import Path

# Add the parent directory to sys.path
parent_dir = str(Path(__file__).resolve().parent.parent)
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

from sklearn.decomposition import PCA



def calc_obsm_pca(adata_file_paths, n_components=50, common_space=False):
    adata_files = {}
    if common_space:
        pca = PCA(n_components=n_components, svd_solver='full')
    for counter, (key, path) in enumerate(adata_file_paths.items()):
        adata = anndata.read_h5ad(path)
        if not common_space:
            pca = PCA(n_components=n_components, svd_solver='full')
            pca.fit(adata.X)
        elif counter == 0:
            pca.fit(adata.X)

        adata.obsm[f'pca_{n_components}'] = pca.transform(adata.X)
        adata_files[key] = adata

    return adata_files


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


def benchmark(adata, latent_adata, batch_key, cell_key, approach):
    nmi_value = nmi(latent_adata, cell_key)
    gc_value = graph_connectivity_score(latent_adata)
    f1_value = isolated_label_f1_score(latent_adata, cell_key, batch_key)
    ari_value = ari_score(adata, cell_key)
    asw_b, asw_c = asw(latent_adata, cell_key, batch_key)
    ebm_value = entropy_batch_mixing(adata, batch_key, n_neighbors=15)
    knn_acc_value = knn_accuracy(latent_adata, cell_key)
    common_metrics = {
        'Approach': approach,
        'NMI': nmi_value,
        'GC': gc_value,
        'ILF1': f1_value,
        'ARI': ari_value,
        'EBM': ebm_value,
        'KNN Acc': knn_acc_value,
        'ASW_B': asw_b,
        'ASW_C': asw_c
    }
    return common_metrics


def tuning_metrics(adata, latent_adata, batch_key, cell_key, epoch, round, approach):
    nmi_value = nmi(latent_adata, cell_key)
    ari_value = ari_score(adata, cell_key)
    ebm_value = entropy_batch_mixing(adata, batch_key, n_neighbors=15)
    asw_b, asw_c = asw(latent_adata, cell_key, batch_key)
    knn_acc_value = knn_accuracy(latent_adata, cell_key)
    common_metrics = {
        'Approach': approach,
        'Epoch': epoch,
        'Round': round,
        'NMI': nmi_value,
        'ARI': ari_value,
        'EBM': ebm_value,
        'KNN Acc': knn_acc_value,
        'ASW_B': asw_b,
        'ASW_C': asw_c
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


def benchmark_all_datasets(fed_data_dir: str, cent_data_dir: str, inclusion: str, n_components: int, batch_key: str,
                           cell_key: str):
    all_metrics = []
    for ds_name in DATASETS:
        if inclusion != "all" and ds_name in ["CellLine", "HumanDendriticCells"]:
            continue
        print(f"Processing {ds_name} {inclusion}...")
        n_clients = 5 if ds_name == "HumanPancreas" else 3 if ds_name == "CellLine" else 2
        fedscgen = anndata.read_h5ad(
            os.path.join(fed_data_dir, ds_name, inclusion, f"BO0-C{n_clients}", "fed_corrected.h5ad"))
        latent_adata = anndata.AnnData(fedscgen.obsm[f'pca_{n_components}'])
        latent_adata.obs = fedscgen.obs
        fedscgen_metrics = benchmark(fedscgen, latent_adata, batch_key, cell_key, "FedscGen")
        fedscgen_metrics["Dataset"] = ds_name
        scgen = anndata.read_h5ad(os.path.join(cent_data_dir, ds_name, inclusion, "corrected.h5ad"))
        latent_adata = anndata.AnnData(scgen.obsm[f'pca_{n_components}'])
        latent_adata.obs = scgen.obs
        scgen_metrics = benchmark(scgen, latent_adata, batch_key, cell_key, "scGen")
        scgen_metrics["Dataset"] = ds_name
        all_metrics.extend([fedscgen_metrics, scgen_metrics])
    df = pd.DataFrame(all_metrics)
    df.to_csv(os.path.join(fed_data_dir, f"fed_cent_metrics-{inclusion}.csv"), sep=",", index=False)


def benchmark_reproduce(fed_data_dir: str, cent_data_dir: str, inclusion: str, n_components: int, batch_key: str,
                        cell_key: str):
    all_metrics = []
    ds_name = "MouseHematopoieticStemProgenitorCells"
    inclusion = "all"
    n_clients = 5 if ds_name == "HumanPancreas" else 3 if ds_name == "CellLine" else 2

    def bench_model(path):
        # adata = anndata.read_h5ad(path)
        adata = calc_obsm_pca({'a': path}, 20)["a"]
        latent_adata = anndata.AnnData(adata.obsm[f'pca_{n_components}'])
        latent_adata.obs = adata.obs
        metrics = benchmark(adata, latent_adata, batch_key, cell_key, "FedscGen")
        return metrics

    fed = bench_model(os.path.join(fed_data_dir, ds_name, inclusion, f"BO0-C2", "fed_corrected.h5ad"))
    # fed = bench_model(os.path.join(fed_data_dir, ds_name, inclusion, f"BO0-C2", f"corrected_10.h5ad"))
    scgen = bench_model(os.path.join(cent_data_dir, ds_name, inclusion, "corrected.h5ad"))
    tun = bench_model(os.path.join(fed_data_dir, ds_name, inclusion, f"BO0-C2", f"corrected_10.h5ad"))
    fed_diff = {k: fed[k] - scgen[k] for k in scgen.keys() if k != "Approach"}
    tun_diff = {k: tun[k] - scgen[k] for k in scgen.keys() if k != "Approach"}
    print(fed_diff)
    print(tun_diff)


def benchmark_batch_out(fed_data_dir, cent_data_dir, n_batches, n_components, batch_key, cell_key):
    """

    Parameters
    ----------
    fed_data_dir
        /home/bba1658/FedscGen/results/scgen/federated/HumanPancreas/all/BO1-C4
    cent_data_dir: str
        /home/bba1658/FedscGen/results/scgen/centralized/HumanPancreas/all
    n_batches
    n_components
    batch_key
    cell_key

    Returns
    -------

    """
    all_metrics = []
    for b in range(n_batches):
        fedscgen = anndata.read_h5ad(os.path.join(fed_data_dir, f"{b}", "fed_corrected_with_new_studies.h5ad"))
        latent_adata = anndata.AnnData(fedscgen.obsm[f'pca_{n_components}'])
        latent_adata.obs = fedscgen.obs
        fedscgen_metrics = benchmark(fedscgen, latent_adata, batch_key, cell_key, "FedscGen")
        fedscgen_metrics["Batch Out"] = b
        scgen = anndata.read_h5ad(os.path.join(cent_data_dir, f"BO{b}", "corrected.h5ad"))
        latent_adata = anndata.AnnData(scgen.obsm[f'pca_{n_components}'])
        latent_adata.obs = scgen.obs
        scgen_metrics = benchmark(scgen, latent_adata, batch_key, cell_key, "scGen")
        scgen_metrics["Batch Out"] = b
        all_metrics.extend([fedscgen_metrics, scgen_metrics])
    df = pd.DataFrame(all_metrics)
    df.to_csv(os.path.join(fed_data_dir, "batchout_metrics.csv"), sep=",", index=False)


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
            files[round] = os.path.join(data_dir, f"E{epoch}", f"corrected_{round}.h5ad")
        files = sorted(files.items(), key=lambda x: x[0])
        h5ad_files[epoch] = files
    all_metrics = []
    sorted_files = sorted(h5ad_files.items(), key=lambda x: x[0])
    for epoch, file in sorted_files:
        for round, path in file:
            adata = anndata.read_h5ad(path)
            latent_adata = anndata.AnnData(adata.obsm[f'pca_{n_components}'])
            latent_adata.obs = adata.obs
            all_metrics.append(
                tuning_metrics(adata, latent_adata, batch_key, cell_key, epoch, round, 'FedscGen'))
    # add scGen
    adata = anndata.read_h5ad(os.path.join(data_dir, "scGen.h5ad"))
    latent_adata = anndata.AnnData(adata.obsm[f'pca_{n_components}'])
    latent_adata.obs = adata.obs
    all_metrics.append(tuning_metrics(adata, latent_adata, batch_key, cell_key, 0, 0, f"scGen"))
    df = pd.DataFrame(all_metrics)
    plot_name = os.path.join(data_dir, "metrics.png")
    df.to_csv(os.path.join(data_dir, "metrics.csv"), sep=",", index=False)
    bar_plot_subplot(df, plot_name)


def load_metrics_and_plot(df_path, plot_name):
    df = pd.read_csv(df_path)
    bar_plot(df, plot_name.replace(".", "-bar."))
    df.drop(columns=["ILS", "ILF1"], inplace=True)
    plot_metrics_with_circles(df, plot_name.replace(".", "-circle."))


def compute_p_values(fedscgen_df, scgen_df):
    """Computes p-values using Wilcoxon signed-rank test for paired samples."""
    p_values = {}

    for metric in fedscgen_df.columns:
        diffs = fedscgen_df[metric] - scgen_df[metric]

        # If all differences are zero, assign p-value of 1 (no significant difference)
        if np.all(diffs == 0):
            p_values[metric] = 1.0
        else:
            try:
                _, p_val = wilcoxon(fedscgen_df[metric], scgen_df[metric], zero_method='wilcox')
                p_values[metric] = p_val
            except ValueError as e:
                print(f"Warning: Wilcoxon test failed for {metric} due to zero differences.")
                p_values[metric] = 1.0  # Assign a non-significant p-value

    return pd.Series(p_values)


def significance_marker(p):
    """Returns significance markers based on p-values."""
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return ""

def calculate_statistical_significance(fed_data_dir, inclusion):
    file_path = os.path.join(fed_data_dir, f"fed_cent_metrics-{inclusion}.csv")
    output_dir = fed_data_dir
    df = pd.read_csv(file_path)
    fedscgen_df = df[df["Approach"] == "FedscGen"].drop(columns=["Approach", "Dataset"])
    scgen_df = df[df["Approach"] == "scGen"].drop(columns=["Approach", "Dataset"])
    p_values = compute_p_values(fedscgen_df, scgen_df)
    significance_df = p_values.apply(significance_marker)
    diff_df = fedscgen_df - scgen_df
    diff_df.to_csv(os.path.join(output_dir, f"performance_diff-{inclusion}.csv"))
    p_values.to_csv(os.path.join(output_dir, f"p_values-{inclusion}.csv"))
    significance_df.to_csv(os.path.join(output_dir, f"significance-{inclusion}.csv"))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate and Plot NMI for a set of adata files.')

    # Add arguments
    parser.add_argument("--data_dir", type=str, help="path to the main data directory",
                        default="/home/bba1658/FedSC/results/scgen/batchout")
    parser.add_argument("--fed_data_dir", type=str, help="path to the main data directory")
    parser.add_argument('--n_components', type=int, default=20, help="Number of components for PCA")
    parser.add_argument("--inclusion", type=str, default="all", choices=["all", "combined", "dropped"])
    parser.add_argument('--cell_key', help='Cell key name.', default="cell_type")
    parser.add_argument('--batch_key', help='Batch key name.', default="batch")
    parser.add_argument("--scenarios", type=str, default="all",
                        choices=["datasets", "batch-out", "snapshots", "tuning", "reproduce", "significance"])
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
            benchmark_all_datasets(args.fed_data_dir, args.data_dir, args.inclusion, args.n_components, args.batch_key,
                                   args.cell_key)
        elif args.scenarios == "batch-out":
            benchmark_batch_out(args.fed_data_dir, args.data_dir, args.n_batches, args.n_components, args.batch_key,
                                args.cell_key)
        elif args.scenarios == "tuning":
            benchmark_tuning(args.data_dir, args.n_components, args.batch_key, args.cell_key)
        elif args.scenarios == "reproduce":
            benchmark_reproduce(args.fed_data_dir, args.data_dir, args.inclusion, args.n_components, args.batch_key,
                                args.cell_key)
        elif args.scenarios == "significance":
            calculate_statistical_significance(args.fed_data_dir, args.inclusion)
        else:
            benchmark_snapshots(args.data_dir, args.n_rounds, args.n_components, args.batch_key, args.cell_key)
