import numpy as np
import pandas as pd
import anndata
import argparse
import os
from scarches.metrics import nmi, entropy_batch_mixing, asw
from pathlib import Path
from sklearn.decomposition import PCA
import sys

# Add the parent directory to sys.path
parent_dir = str(Path(__file__).resolve().parent.parent)
if parent_dir not in sys.path:
    sys.path.append(parent_dir)
from utils import (graph_connectivity_score, isolated_label_f1_score, ari_score, bar_plot, plot_metrics_with_circles,
                   compute_ils, knn_accuracy, DATASETS, bar_plot_subplot)


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


def benchmark(adata, latent_adata, batch_key, cell_key):
    nmi_value = nmi(latent_adata, cell_key)
    gc_value = graph_connectivity_score(latent_adata)
    f1_value = isolated_label_f1_score(latent_adata, cell_key, batch_key)
    ari_value = ari_score(adata, cell_key)
    asw_b, asw_c = asw(latent_adata, cell_key, batch_key)
    ebm_value = entropy_batch_mixing(adata, batch_key, n_neighbors=15)
    knn_acc_value = knn_accuracy(latent_adata, cell_key)
    common_metrics = {
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


def benchmark_dataset(file_path, n_components, batch_key, cell_key):
    results = anndata.read_h5ad(file_path)
    latent_adata = anndata.AnnData(results.obsm[f'pca_{n_components}'])
    latent_adata.obs = results.obs
    metrics = benchmark(results, latent_adata, batch_key, cell_key)
    return metrics



def find_all_corrected_files_path(res_dir, approach):
    inclusion_scenarios = ["all", "combined", "dropped"]
    df = pd.DataFrame(columns=["Seed",
                               "Epoch",
                               "Round",
                               "File",
                               "Inclusion",
                               "BatchOut",
                                 "Batch",
                                    "N_Clients"
                               ]).astype({
                                    "Seed": int,
                                    "Epoch": int,
                                    "Round": int,
                                    "File": str,
                                    "Inclusion": str,
                                    "BatchOut": int,
                                    "Batch": int,
                                    "N_Clients": int
                                })
    # find all "corrected.h5ad"
    corrected_files = [file for file in Path(res_dir).rglob("*.h5ad") if file.is_file()]
    for file in corrected_files:
        inclusion = "all"
        for scenario in inclusion_scenarios:
            if scenario in str(file):
                inclusion = scenario
                break
        basename = Path(file).stem
        seed = int(str(file).split('seed_')[1].split('/')[0]) if 'seed_' in str(file) else 42
        epoch, round, batch_out, batch, n_clients = np.nan, np.nan, np.nan, np.nan, np.nan
        parent_dir = Path(file).parent
        if "corrected_" in basename and 'with_new_studies' not in basename:
            try:
                round = int(basename.split("corrected_")[1].split(".h5ad")[0])
            except:
                raise ValueError(f"Error in finding rounds in file: {file}")
            try:
                epoch = int(parent_dir.name.split("E")[1])
            except:
                raise ValueError(f"Error in finding epochs in file: {file}")
        if '/BO' in str(file):
            try:
                if approach == "scgen":
                    batch_out = 1
                    batch = int(str(file).split('BO')[1][0])
                else:
                    batch_out = int(str(file).split('BO')[1][0])

                    n_clients = str(file).split(f'BO{batch_out}-C')[1][0]
                    if parent_dir.name.isdigit():
                        batch = int(parent_dir.name)
            except:
                raise ValueError(f"Error in finding batch out in file: {file}")
        df = pd.concat([df, pd.DataFrame({"Seed": seed,
                                            "Epoch": epoch,
                                            "Round": round,
                                            "File": str(file),
                                            "Inclusion": inclusion,
                                            "BatchOut": batch_out,
                                            "Batch": batch,
                                            "N_Clients": n_clients
                                          },
                                         index=[0])])
    return df

def benchmark_all(data_dir: str, approach: str, n_components, batch_key, cell_key, tuning=False):
    approaches = {"scgen": "scGen", "fedscgen": "FedscGen", "fedscgen-smpc": "FedscGen-SMPC"}
    if tuning:
        assert approach == "fedscgen", "Only FedscGen approach is supported for tuning"
        res_dir = data_dir
        approach = "FedscGen"
    else:
        assert approach in approaches.keys(), "Approach must be one of {}".format(approaches.keys())
        res_dir = os.path.join(data_dir, approach)
        approach = approaches[approach]
    output_file = os.path.join(res_dir, f"benchmark_metrics.csv")
    file_exists = os.path.exists(output_file) and os.stat(output_file).st_size > 0
    if file_exists:
        results_df = pd.read_csv(output_file)
    else:
        results_df = pd.DataFrame(columns=["Seed",
                                           "Epoch",
                                           "Round",
                                           "File",
                                           "Inclusion",
                                           "Dataset",
                                           "Approach",
                                           "BatchOut",
                                             "Batch",
                                                "N_Clients",
                                             "NMI",
                                                "GC",
                                                "ILF1",
                                           "ARI",
                                             "EBM",
                                                "KNN Acc",
                                                "ASW_B",
                                                "ASW_C"
                                           ])
    for ds_name in DATASETS:
        if ds_name in os.listdir(res_dir):
            df = find_all_corrected_files_path(os.path.join(res_dir, ds_name))
            df["Dataset"] = ds_name
            df["Approach"] = approach
            for index, row in df.iterrows():
                if row['File'] in results_df["File"].values:
                    continue
                try:
                    print(f"[BENCHMARK] ==> {row['File']}")
                    # metric = benchmark_dataset(row['File'], n_components, batch_key, cell_key)
                    metric = {"NMI": 0.0, "GC": 0.0, "ILF1": 0.0, "ARI": 0.0, "EBM": 0.0, "KNN Acc": 0.0, "ASW_B": 0.0,}
                    new_row_df = pd.DataFrame([{**row.to_dict(), **metric}])
                    write_header = not os.path.exists(output_file) or os.stat(output_file).st_size == 0
                    new_row_df.to_csv(output_file, mode="a", sep=",", index=False, header=write_header)
                except Exception as e:
                    print(f"Error processing dataset: {ds_name}, file: {row['File']}, seed: {row['Seed']}. Error: {e}")
                    continue


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
    parser.add_argument("--fed_data_dir", type=str, help="path to the main data directory")
    parser.add_argument('--n_components', type=int, default=20, help="Number of components for PCA")
    parser.add_argument("--inclusion", type=str, default="all", choices=["all", "combined", "dropped"])
    parser.add_argument('--cell_key', help='Cell key name.', default="cell_type")
    parser.add_argument('--batch_key', help='Batch key name.', default="batch")
    parser.add_argument("--scenarios", type=str, default="all",
                        choices=["approach", "tuning", "reproduce"])
    parser.add_argument("--approach", type=str, default="fedscgen", choices=["fedscgen", "scgen", "fedscgen-smpc"])
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
        if args.scenarios == "approach":
            benchmark_all(args.data_dir, args.approach, args.n_components, args.batch_key,
                          args.cell_key)
        elif args.scenarios == "tuning":
            benchmark_all(args.data_dir, "fedscgen", args.n_components, args.batch_key, args.cell_key, tuning=True)
        elif args.scenarios == "reproduce":
            benchmark_reproduce(args.fed_data_dir, args.data_dir, args.inclusion, args.n_components, args.batch_key,
                                args.cell_key)
        else:
            benchmark_snapshots(args.data_dir, args.n_rounds, args.n_components, args.batch_key, args.cell_key)
