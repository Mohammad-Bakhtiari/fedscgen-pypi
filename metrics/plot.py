import copy
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.colors as colors
import pandas as pd
import numpy as np
import os
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import itertools
import argparse
from utils import DATASETS_COLORS, DATASETS_MARKERS, DATASETS_ACRONYM, DATASETS
import glob
import sys
from pathlib import Path
from scipy import stats
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Add the parent directory to sys.path
parent_dir = str(Path(__file__).resolve().parent.parent)
if parent_dir not in sys.path:
    sys.path.append(parent_dir)



def read_metrics_files(data_dir, filename="metrics.csv"):
    """ find all files named as `metrics.csv` in the data_dir and read them into a single dataframe

    Parameters
    ----------
    data_dir

    Returns
    -------

    """
    all_metrics = []
    for dataset, acr in zip(DATASETS, DATASETS_ACRONYM):
        metrics_file = os.path.join(data_dir, dataset, filename)
        if os.path.exists(metrics_file):
            metrics = pd.read_csv(metrics_file)
            metrics["Dataset"] = acr
            all_metrics.append(metrics)
    return pd.concat(all_metrics)


def plot_metrics_with_heatmaps(df, metric_keys, plot_name):
    scGen = df[df["Approach"] == "scGen"]
    df = df[df["Approach"] == "FedscGen"]
    df.drop(columns=["Approach"], inplace=True)

    dataset_keys = df['Dataset'].unique().tolist()

    fig, axs = plt.subplots(len(metric_keys), len(dataset_keys), figsize=(len(dataset_keys) * 4, df.Round.max() * 2),
                            squeeze=False)
    cmap = cm.get_cmap('viridis')
    norm = colors.Normalize(vmin=df[metric_keys].min().min(), vmax=df[metric_keys].max().max())
    mappable = cm.ScalarMappable(norm=norm, cmap=cmap)

    for j, dataset in enumerate(dataset_keys):
        for i, metric in enumerate(metric_keys):
            ax = axs[i][j]
            data = df[df['Dataset'] == dataset]
            # pivot = data.pivot(index='Round', columns='Epoch', values=metric)
            pivot = data.pivot(index='Epoch', columns='Round', values=metric)
            sns.heatmap(pivot, ax=ax, cmap='viridis', cbar=False, vmin=0, vmax=1)

            scgen_color_value = scGen[scGen['Dataset'] == dataset][metric].values[0]
            scgen_color = mappable.to_rgba(scgen_color_value)

            # Create a legend for each subplot
            legend_elements = [plt.Line2D([0], [0], marker='o', color='w', label='scGen',
                                          markerfacecolor=scgen_color, markersize=15, markeredgewidth=1,
                                          markeredgecolor='black')]
            ax.legend(handles=legend_elements, loc=(0.6, 1.02))
            ax.set_xticklabels(ax.get_xticklabels(), fontsize=14)
            ax.set_yticklabels(ax.get_yticklabels(), fontsize=14)
            ax.grid(False)
            if j == 0:  # Optionally, label rows with dataset names
                ax.text(-0.5, 0.38, metric, transform=ax.transAxes, fontsize=24, va='bottom', ha='left', rotation=90)
                ax.set_ylabel("Epochs", fontsize=18)
            else:
                ax.set_ylabel('')
                # set off the y-axis
                ax.yaxis.set_visible(False)
            if i == 0:  # Optionally, label columns with metric names
                ax.set_title(dataset, fontsize=18, loc='left')
            if i < len(metric_keys) - 1:
                ax.xaxis.set_visible(False)
            else:
                ax.set_xlabel("Rounds", fontsize=18)

    # Adjust layout for colorbar and legend
    plt.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.805, 0.28, 0.01, 0.5])
    cbar = plt.colorbar(mappable, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=14)

    # Adjust layout to make space for legends
    plt.subplots_adjust(right=0.8)
    plt.savefig(plot_name, dpi=300)


def plot_tuning_heatmap(dataset_keys, df, metric_keys, plot_name, scGen):
    fig, axs = plt.subplots(len(metric_keys), len(dataset_keys),
                            figsize=(len(dataset_keys) * 4, len(df['Epoch'].unique()) * 2),
                            squeeze=True)
    fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.01, hspace=0.01)
    # Ensure symmetric vmin and vmax for color normalization around zero
    max_abs_value = max(abs(df[metric_keys].min().min()), abs(df[metric_keys].max().max()))
    max_abs_value = 1
    norm = colors.Normalize(vmin=-max_abs_value, vmax=max_abs_value)
    cmap = cm.get_cmap('RdBu')
    mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
    for j, dataset in enumerate(dataset_keys):
        for i, metric in enumerate(metric_keys):
            ax = axs[i][j]
            data = df[df['Dataset'] == dataset]
            scgen_ds_metric = scGen[scGen["Dataset"] == dataset][metric].values[0]
            # Find the index of the maximum value
            max_index = data[metric].idxmax()
            best_epoch = data.loc[max_index, "Epoch"]
            best_round = data.loc[max_index, "Round"]
            data[metric] = data[metric].apply(lambda x: x - scgen_ds_metric)
            max_value_diff = data.loc[max_index, metric]

            pivot = data.pivot(index='Epoch', columns='Round', values=metric)
            # abs_max = max(abs(pivot.min().min()), abs(pivot.max().max()))
            abs_max = 1
            cmap = 'RdBu'
            sns.heatmap(pivot, ax=ax, cmap=cmap, cbar=False, center=0, vmin=-abs_max, vmax=abs_max, square=True)
            ax.text(best_round - 0.5, best_epoch - 0.5, f"{max_value_diff:.2f}", ha='center', va='center', fontsize=12)

            ax.set_xticklabels(ax.get_xticklabels(), fontsize=14)
            ax.set_yticklabels(ax.get_yticklabels(), fontsize=14)
            ax.grid(False)
            if j == 0:
                ax.text(-0.5, 0.38, metric, transform=ax.transAxes, fontsize=24, va='bottom', ha='left',
                        rotation=90)
                ax.set_ylabel("Epochs", fontsize=18)
            else:
                ax.set_ylabel('')
                ax.yaxis.set_visible(False)
            if i == 0:
                ax.set_title(dataset, fontsize=18, loc='left')
            if i < len(metric_keys) - 1:
                ax.xaxis.set_visible(False)
            else:
                ax.set_xlabel("Rounds", fontsize=18)
    # Adjust layout for colorbar and legend
    plt.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.805, 0.28, 0.01, 0.5])
    cbar = plt.colorbar(mappable, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=14)
    plt.savefig(plot_name, dpi=300)


def plot_bo_hitmap(df, plt_name, dpi, font_size=20, tick_size=14, cell_size=1, cbar_font_size=14, tick_label_size=26):
    # Determine the figure size based on the number of columns and rows
    num_cols, num_rows = df.shape
    fig_width = cell_size * num_rows
    fig_height = cell_size * num_cols

    # Plotting the heatmap
    fig = plt.figure(figsize=(fig_width + 3, fig_height + 3))
    fig.subplots_adjust(left=0.05, right=0.8, top=0.95, bottom=0.2, wspace=0.01, hspace=0.01)
    plt.grid(False)
    # abs_max = max(abs(df.values.min()), abs(df.values.max()))
    abs_max = 0.5
    ax = sns.heatmap(df, annot=True, cmap='RdBu', vmin=-abs_max, vmax=abs_max, center=0, annot_kws={"size": 26},
                     square=True, cbar=False)
    plt.xticks(fontsize=tick_size)
    ax.yaxis.set_visible(False)
    ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=30, rotation=45, ha='right')
    norm = colors.Normalize(vmin=-abs_max, vmax=abs_max)
    mappable = plt.cm.ScalarMappable(cmap='RdBu', norm=norm)
    cbar_ax = fig.add_axes([0.81, 0.27, 0.02, 0.6])
    cbar = plt.colorbar(mappable, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=26)
    plt.savefig(plt_name, dpi=dpi)


def find_best_round_epoch(dataset_keys, df, metric_keys, scGen):
    dfs = []
    for dataset in dataset_keys:
        data = df[df['Dataset'] == dataset]
        for metric in metric_keys:
            scgen_ds_metric = scGen[scGen["Dataset"] == dataset][metric].values[0]
            data[metric] = data[metric].apply(lambda x: x - scgen_ds_metric)
        dfs.append(data)
    df = pd.concat(dfs)
    res_min = {}
    for r in df['Round'].unique():
        for e in df['Epoch'].unique():
            data = df[(df['Round'] == r) & (df['Epoch'] == e)]
            res_min[f"E{e}R{r}"] = data[metric_keys].min().min()
    res_max = sorted(res_min.items(), key=lambda x: x[1], reverse=True)
    print(res_max[:10])
    count_max = {}
    for k, v in res_max:
        e = int(k[1:].split("R")[0])
        r = int(k.split("R")[1])
        temp = df[(df['Round'] == r) & (df['Epoch'] == e)]
        count_max[k] = (temp[metric_keys] >= 0).sum().sum()
    count_max = sorted(count_max.items(), key=lambda x: x[1], reverse=True)
    print(count_max[:10])
    print(len(count_max))
    er = df[(df["Epoch"] == 2) & (df["Round"] == 8)]
    print((er[metric_keys] >= 0).sum())
    print((er[metric_keys] < 0).sum())


def read_kbet_score(scores_path):
    df = pd.read_csv(scores_path)
    df["stats"] = df['Unnamed: 0']
    df['acceptance_rate'] = 1 - df['kBET.observed']
    df = df.groupby(["approach", "k_value"])["acceptance_rate"].median().reset_index()
    return df


def acceptance_plot(df, plot_name):
    # k_values = df['k_values'].unique()

    for approach, color in zip(["scGen", "FedscGen"], ['blue', 'red']):
        subset = df[df["approach"] == approach]
        plt.plot(subset['k_value'], subset['acceptance_rate'], label=approach, marker='o', color=color)

    plt.xlabel('k-value')
    plt.ylabel('Acceptance Rate')
    plt.xticks(df.k_value.unique())
    plt.legend()
    plt.title('kBET Acceptance Rates across k-values')
    plt.savefig(plot_name)
    plt.close()


def read_tuning_res(data_dir, read_all=False):
    if read_all:
        df = read_metrics_files(data_dir, filename="metrics.csv")
        df.to_csv(os.path.join(data_dir, "all_metrics.csv"), index=False)
    df = pd.read_csv(os.path.join(data_dir, "all_metrics.csv"))
    metric_keys = ["ARI", "NMI", "EBM", "ASW_B", "ASW_C", "KNN Acc"]
    plot_name = os.path.join(data_dir, "tuning-diff.png")
    scGen = df[df["Approach"] == "scGen"]
    df = df[df["Approach"] == "FedscGen"]
    df.drop(columns=["Approach"], inplace=True)
    dataset_keys = df['Dataset'].unique().tolist()
    find_best_round_epoch(dataset_keys, copy.deepcopy(df), metric_keys, scGen)
    plot_tuning_heatmap(dataset_keys, df, metric_keys, plot_name, scGen)


def acceptance_diff_plot(df, plot_name):
    unique_datasets = df["Dataset"].unique()

    # Group by dataset, approach, and k_value to compute the median acceptance rate
    grouped_df = df.groupby(["Dataset", "approach", "k_value"])["acceptance_rate"].median().reset_index()

    # Pivot the dataframe to get acceptance rates for each approach in separate columns
    pivoted_df = grouped_df.pivot_table(index=["Dataset", "k_value"], columns="approach",
                                        values="acceptance_rate").reset_index()

    # Calculate the difference in acceptance rates
    pivoted_df["difference"] = pivoted_df["FedscGen"] - pivoted_df["scGen"]
    df = pivoted_df

    # Create a Seaborn line plot
    plt.figure(figsize=(6, 4))
    for i, dataset in enumerate(unique_datasets):
        subset = df[df["Dataset"] == dataset]
        c_num = i + 1 if i < 3 else i + 1
        sns.scatterplot(data=subset, x='k_value', y='difference', label=dataset,
                        color=DATASETS_COLORS[dataset],
                        marker=DATASETS_MARKERS[dataset], s=100)

    # Adjust legend position to be outside the plot
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=20, prop={"size": 12})
    plt.axhline(y=0, color='gray', linestyle='-')
    plt.xlabel('k-value', fontsize=16)
    plt.ylabel('Acceptance Rate Difference', fontsize=16)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.xticks([0.05, 0.1, 0.15, 0.2, 0.25], fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.savefig(plot_name, dpi=1000)
    plt.close()


def read_kbet_file(hp_kbet_file):
    kbet_df = pd.read_csv(hp_kbet_file)
    kbet_df.drop(columns=["Unnamed: 0", "filename"], inplace=True)
    kbet_df["approach"] = kbet_df["approach"].apply(lambda x: "scGen" if x == "scgen" else "FedscGen")
    kbet_df['acceptance_rate'] = 1 - kbet_df['kBET.observed']
    grouped_df = kbet_df.groupby(["approach", "k_value"])["acceptance_rate"].median().reset_index()
    pivoted_df = grouped_df.pivot_table(index=["k_value"], columns="approach",
                                        values="acceptance_rate").reset_index()
    pivoted_df["difference"] = pivoted_df["FedscGen"] - pivoted_df["scGen"]
    kbet_diff = pivoted_df["difference"].mean().round(2)
    return kbet_diff


def read_batchout(bo_metrics_file, hp_kbet_file, bo_kbet_dir, all_ds_metrics_file, output_dir):
    hp = get_hp_metrics_kbet_diff(all_ds_metrics_file, hp_kbet_file)
    difference_df = get_bo_metrics_kbet_diff(bo_kbet_dir, bo_metrics_file)
    difference_df = pd.concat([difference_df, hp])
    difference_df.to_csv(f"{output_dir}/bo-metrics.csv", index=True)
    difference_df = pd.read_csv(f"{output_dir}/bo-metrics.csv")
    difference_df.set_index("Batch Out", inplace=True)
    plot_bo_hitmap(difference_df, f"{output_dir}/bo-hitmap.png", dpi=300)


def get_bo_metrics_kbet_diff(bo_kbet_dir, bo_metrics_file):
    kbet_diff_values = {}
    for bo in range(5):
        file_path = os.path.join(bo_kbet_dir, str(bo), 'kBET_summary_results.csv')
        kbet_diff_values[bo + 1] = read_kbet_file(file_path)
    bo_metrics_df = pd.read_csv(bo_metrics_file)
    bo_metrics_df["Batch Out"] = bo_metrics_df["Batch Out"].apply(lambda x: int(x + 1))
    if "ILS" in bo_metrics_df.columns:
        bo_metrics_df.drop(columns=["ILS"], inplace=True)
    scgen_df = bo_metrics_df[bo_metrics_df['Approach'] == 'scGen'].set_index('Batch Out')
    scgen_df.drop(columns=["Approach"], inplace=True)
    fedscgen_df = bo_metrics_df[bo_metrics_df['Approach'] == 'FedscGen'].set_index('Batch Out')
    fedscgen_df.drop(columns=["Approach"], inplace=True)
    difference_df = fedscgen_df.subtract(scgen_df).round(2)
    difference_df["kBET"] = kbet_diff_values.values()
    return difference_df


def get_hp_metrics_kbet_diff(all_ds_metrics_file, hp_kbet_file):
    all_ds_metrics_df = pd.read_csv(all_ds_metrics_file)
    hp = all_ds_metrics_df[all_ds_metrics_df["Dataset"] == "HumanPancreas"]
    hp.drop(columns=["Dataset"], inplace=True)
    fedscgen_hp = hp[hp["Approach"] == "FedscGen"].reset_index()
    fedscgen_hp.drop(columns=["Approach", "index"], inplace=True)
    scgen_hp = hp[hp["Approach"] == "scGen"].reset_index()
    scgen_hp.drop(columns=["Approach", "index"], inplace=True)
    hp = fedscgen_hp.subtract(scgen_hp).round(2)
    kbet_value_diff = read_kbet_file(hp_kbet_file)
    hp["kBET"] = kbet_value_diff
    hp["Batch Out"] = None
    hp.set_index("Batch Out", inplace=True)
    return hp


def read_kbet(data_dir):
    for inclusion in ["all", "combined", "dropped"]:
        print("Inclusion:", inclusion)
        all_kbets = []
        for dataset, acr in zip(DATASETS, DATASETS_ACRONYM):
            if dataset == "MouseBrain" or (inclusion != "all" and dataset in ["CellLine", "HumanDendriticCells"]):
                continue
            n_clients = 5 if dataset == "HumanPancreas" else 3 if dataset == "CellLine" else 2
            kbet_file = os.path.join(data_dir, dataset, inclusion, f"BO0-C{n_clients}", "kBET_summary_results.csv")
            if os.path.exists(kbet_file):
                kbet_df = pd.read_csv(kbet_file)
                kbet_df.drop(columns=["Unnamed: 0", "filename"], inplace=True)
                kbet_df["Dataset"] = acr
                kbet_df["approach"] = kbet_df["approach"].apply(lambda x: "scGen" if x == "scgen" else "FedscGen")
                all_kbets.append(kbet_df)
            else:
                raise FileNotFoundError(f"{kbet_file} is not exist")
        df = pd.concat(all_kbets)
        df.sort_values(by="Dataset", inplace=True)
        df['acceptance_rate'] = 1 - df['kBET.observed']

        acceptance_diff_plot(df, os.path.join(data_dir, f'kBET_Acceptance_Rates-{inclusion}.png'))


def read_scenarios_metrics(data_dir):
    all_scenarios = []
    for inclusion in ["all", "combined", "dropped"]:
        df = get_scenario_metrics_diff(data_dir, inclusion, skip_datasets=["CellLine", "HumanDendriticCells"])
        all_scenarios.append(df)
    df = pd.concat(all_scenarios, ignore_index=True)
    # df = pd.read_csv(os.path.join(data_dir, "datasets-metrics-scenarios.csv"))
    df["inclusion"] = df["inclusion"].apply(lambda x: x.capitalize())
    # subplot heatmap of each dataset for all metrics vs inclusion
    fig, axs = plt.subplots(2, 3, figsize=(30, 8), squeeze=True)
    fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.01, hspace=0.01)
    i = 0
    j = 0
    abs_max = max([abs(df.drop(columns=["Dataset", "inclusion"]).min().min()),
                   abs(df.drop(columns=["Dataset", "inclusion"]).max().max())]
                  )
    abs_max = 1
    norm = colors.Normalize(vmin=-abs_max, vmax=abs_max)
    cmap = cm.get_cmap('RdBu')
    mappable = cm.ScalarMappable(norm=norm, cmap=cmap)

    for dataset in ["HP", "MB", "MR", "MCA", "PBMC", "MHSPC"]:
        subset = df[df["Dataset"] == dataset].set_index("inclusion")
        subset.drop(columns=["Dataset"], inplace=True)
        if i == 3:
            i = 0
            j += 1
        ax = axs[j][i]

        sns.heatmap(subset, annot=True, cmap='RdBu', center=0, ax=ax, cbar=False, annot_kws={"size": 16}, square=True,
                    vmin=-abs_max, vmax=abs_max)
        ax.set_title(dataset, fontsize=30)
        ax.grid(False)
        if i == 0:
            ax.set_ylabel("Inclusion", fontsize=30)
            ax.set_yticklabels(ax.get_yticklabels(), fontsize=20, rotation=0, ha='right')
        else:
            ax.set_ylabel('')
            ax.yaxis.set_visible(False)
        if j == 1:
            ax.set_xlabel("Metrics", fontsize=30)
            ax.set_xticklabels(ax.get_xticklabels(), fontsize=20, rotation=45, ha='right')
        else:
            ax.set_xlabel('')
            ax.xaxis.set_visible(False)
        i += 1
    # plot cbar
    plt.subplots_adjust(wspace=0.05, hspace=0.1)
    cbar_ax = fig.add_axes([0.91, 0.25, 0.015, 0.5])
    cbar = plt.colorbar(mappable, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=20)
    plt.savefig(os.path.join(data_dir, "datasets-metrics-scenarios.png"), dpi=300)


def read_datasets_metrics(data_dir):
    df = get_scenario_metrics_diff(data_dir, inclusion="all")
    df.to_csv(os.path.join(data_dir, "datasets-metrics-all.csv"), index=False)
    df = pd.read_csv(os.path.join(data_dir, "datasets-metrics-all.csv"))
    df.set_index("Dataset", inplace=True)
    df.drop(columns=["inclusion"], inplace=True)
    fig, axs = plt.subplots(1, 1, figsize=(16, 14), squeeze=True)
    fig.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.1, wspace=0.01, hspace=0.01)
    abs_max = max([abs(df.min().min()),
                   abs(df.max().max())]
                  )
    abs_max = 1
    sns.heatmap(df, annot=True, cmap='RdBu', center=0, ax=axs, cbar=False, annot_kws={"size": 28}, square=True,
                vmin=-abs_max, vmax=abs_max)
    axs.set_title("Performance difference: FedscGen - scGen", fontsize=36)
    axs.grid(False)
    axs.set_ylabel("Datasets", fontsize=30)
    axs.set_yticklabels(axs.get_yticklabels(), fontsize=30, rotation=0, ha='right')
    axs.set_xlabel("Metrics", fontsize=30)
    axs.set_xticklabels(axs.get_xticklabels(), fontsize=30, rotation=45, ha='right')
    norm = colors.Normalize(vmin=-abs_max, vmax=abs_max)
    mappable = plt.cm.ScalarMappable(cmap='RdBu', norm=norm)
    cbar_ax = fig.add_axes([0.91, 0.25, 0.02, 0.5])
    cbar = plt.colorbar(mappable, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=24)
    plt.savefig(os.path.join(data_dir, "datasets-metrics-all.png"), dpi=300)


def get_scenario_metrics_diff(data_dir, inclusion, skip_datasets=None):
    kbet_df = {}
    for dataset, acr in zip(DATASETS, DATASETS_ACRONYM):
        if inclusion != "all" and dataset in ["CellLine", "HumanDendriticCells"]:
            continue
        if skip_datasets:
            if dataset in skip_datasets:
                continue
        n_clients = 5 if dataset == "HumanPancreas" else 3 if dataset == "CellLine" else 2
        kbet_file = os.path.join(data_dir, dataset, inclusion, f"BO0-C{n_clients}", "kBET_summary_results.csv")
        kbet_df[acr] = read_kbet_file(kbet_file)
    kbet_df = pd.DataFrame(data=kbet_df.values(), columns=["kBET"], index=kbet_df.keys())
    metrics_df = pd.read_csv(os.path.join(data_dir, f"fed_cent_metrics-{inclusion}.csv"))
    if inclusion != "all":
        metrics_df = metrics_df[~metrics_df["Dataset"].isin(["CellLine", "HumanDendriticCells"])]
    metrics_df["Dataset"] = metrics_df["Dataset"].apply(lambda x: DATASETS_ACRONYM[DATASETS.index(x)])
    fedscgen = metrics_df[metrics_df["Approach"] == "FedscGen"].set_index("Dataset")
    fedscgen.drop(columns=["Approach"], inplace=True)
    scgen = metrics_df[metrics_df["Approach"] == "scGen"].set_index("Dataset")
    scgen.drop(columns=["Approach"], inplace=True)
    metric = fedscgen.subtract(scgen).round(2)
    metric["kBET"] = kbet_df
    metric["inclusion"] = inclusion
    metric.reset_index(inplace=True)
    return metric


def get_classification_accuracy(data_dir):
    dfs = []
    for filename in glob.glob(data_dir + "/classification_acc_*.csv"):
        df = pd.read_csv(filename, index_col=None, header=0)
        mean_df = df.groupby('Epoch').mean().reset_index()
        dfs.append(mean_df)
    if len(dfs) == 0:
        return None
    result = pd.concat(dfs, axis=0, ignore_index=True)
    dir_mean_acc_by_epoch = result.groupby('Epoch')['Accuracy'].mean()
    dir_mean_auc_by_epoch = result.groupby('Epoch')['AUC'].mean()

    # Get highest accuracy and its corresponding epoch
    max_acc_epoch = dir_mean_acc_by_epoch.idxmax()
    max_acc = dir_mean_acc_by_epoch.max()

    # Get highest AUC and its corresponding epoch
    max_auc_epoch = dir_mean_auc_by_epoch.idxmax()
    max_auc = dir_mean_auc_by_epoch.max()
    return max_acc_epoch, max_acc, max_auc_epoch, max_auc


def read_classification(data_dir):
    columns = ['Dataset', 'Inclusion', 'Model', "Batch Out", 'Approach', 'Max Accuracy', 'Max Accuracy Epoch',
               'Max AUC', 'Max AUC Epoch']
    results_df = pd.DataFrame(columns=columns)

    for dataset, acr in zip(DATASETS, DATASETS_ACRONYM):
        for model in ["mlp-norm", "knn"]:
            for approach in ["scGen", "FedscGen"]:
                if approach == "scGen":
                    r_dir = os.path.join(data_dir, "centralized", dataset, "all", "classification", model)
                else:
                    n_clients = 5 if dataset == "HumanPancreas" else 3 if dataset == "CellLine" else 2
                    r_dir = os.path.join(data_dir, "federated", dataset, "all", f"BO0-C{n_clients}", "classification",
                                         model)

                metrics = get_classification_accuracy(r_dir)
                if metrics:
                    max_acc_epoch, max_acc, max_auc_epoch, max_auc = metrics
                else:
                    print(f"Metrics file not found in {r_dir}")
                    continue

                row = pd.DataFrame([{
                    'Dataset': acr,
                    'Inclusion': "all",
                    "Model": "MLP" if model == "mlp-norm" else model.upper(),
                    'Approach': approach,
                    'Max Accuracy': max_acc,
                    'Max Accuracy Epoch': max_acc_epoch,
                    'Max AUC': max_auc,
                    'Max AUC Epoch': max_auc_epoch
                }]
                )
                results_df = pd.concat([results_df, row], ignore_index=True)
    print(results_df.head(len(results_df)))
    results_df.to_csv(os.path.join(data_dir, "latent_acc_diff.csv"))

    # Calculate the difference in accuracy between FedscGen and scGen
    diff_df = results_df.pivot_table(index=['Dataset', 'Model'], columns='Approach',
                                     values='Max Accuracy').reset_index()

    diff_df['Accuracy Difference'] = diff_df['FedscGen'] - diff_df['scGen']
    print(diff_df)
    # Plotting
    sns.set(style="white", context="talk")
    fig, ax = plt.subplots(figsize=(6, 6))

    markers = {'MLP': 'o', 'KNN': '*'}
    # use palette='viridis' to color the markers
    palette = sns.color_palette('viridis', len(diff_df['Dataset'].unique()))

    for i, model in enumerate(diff_df["Model"].unique()):
        model_df = diff_df[diff_df["Model"] == model]
        ax.scatter(model_df["Dataset"], model_df["Accuracy Difference"], marker=markers[model], label=model,
                   facecolors='none', edgecolors=palette[i])

    ax.set_xlabel('Dataset')
    ax.set_ylabel('Accuracy Difference')
    ax.axhline(0, color='gray', linestyle='--')
    ax.tick_params(axis='x', rotation=45)
    ax.legend()

    # Format y-axis numbers to two decimal places
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.2f}'))
    ax.set_ylim(-0.15, 0.15)
    # set ytick labels with step of 0.05
    ax.set_yticks(np.arange(-0.15, 0.16, 0.05))

    plt.tight_layout()
    plt.savefig(os.path.join(data_dir, "latent_accuracy_comparison.png"), dpi=300)


def set_fontsize(ax, y_label, font_size, tick_fontsize):
    ax.set_xlabel('', fontsize=font_size)
    ax.set_ylabel(y_label, fontsize=font_size)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=tick_fontsize)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=tick_fontsize)


def read_lisi(data_dir):
    inclusions = ["all", "dropped", "combined"]
    color_palette = sns.color_palette("viridis", 3)
    for inclusion in inclusions:

        # List to hold data from all datasets under a particular target folder.
        aggregated_data = []

        # Define hatch patterns for file_name
        hatch_patterns = ['/', '\\', '|']

        # Loop through each dataset and read its CSV.
        for ind, dataset in enumerate(DATASETS):
            if inclusion != "all" and dataset in ["CellLine", "HumanDendriticCells"]:
                continue
            n_clients = 5 if dataset == "HumanPancreas" else 3 if dataset == "CellLine" else 2
            csv_path = os.path.join(data_dir, dataset, inclusion, f"BO0-C{n_clients}", 'lisi_results.csv')
            if os.path.exists(csv_path):
                df = pd.read_csv(csv_path)
                df.rename(columns={"file_name": "Approach"}, inplace=True)
                df["Approach"] = df["Approach"].apply(
                    lambda x: "FedscGen" if x == "fedscgen" else "scGen" if x == "scgen" else "Raw")
                df["Dataset"] = DATASETS_ACRONYM[ind]
                aggregated_data.append(df)
            else:
                raise FileNotFoundError(f"{csv_path} is not exist")

        # Concatenate all dataframes to get a single dataframe for the current target.
        aggregated_df = pd.concat(aggregated_data, ignore_index=True)
        # order the dataframe raw, scgen, fedscgen
        aggregated_df['Approach'] = pd.Categorical(aggregated_df['Approach'], categories=["Raw", "scGen", "FedscGen"],
                                                   ordered=True)
        # aggregated_df.sort_values(by="Dataset", inplace=True)
        plt.figure(figsize=(18, 6))
        approaches = aggregated_df["Approach"].unique()
        file_name_to_color = {file_name: color for file_name, color in zip(approaches, color_palette)}
        plt.subplot(1, 2, 1)
        ax1 = sns.boxplot(x='Dataset', y='batch', hue='Approach', data=aggregated_df,
                          palette=file_name_to_color)
        ax1.legend(fontsize=16)

        # Boxplot for 'cell_type'
        plt.subplot(1, 2, 2)
        ax2 = sns.boxplot(x='Dataset', y='cell_type', hue='Approach', data=aggregated_df,
                          palette=file_name_to_color)
        ax2.legend_.remove()
        set_fontsize(ax1, 'iLISI', font_size=22, tick_fontsize=20)
        set_fontsize(ax2, 'cLISI', font_size=22, tick_fontsize=20)
        # Add hatch patterns for models (file_name)
        for i, patches in enumerate(zip(ax1.artists, ax2.artists)):
            hatch = hatch_patterns[i % len(hatch_patterns)]
            patches[0].set_hatch(hatch)
            patches[1].set_hatch(hatch)

        # Adjust layout and save
        plt.tight_layout()
        plt.savefig(f'{data_dir}/lisi_{inclusion}.png', dpi=1000)
        plt.close()

def plot_accuracy_diff(df, plot_dir):
    # Melt the DataFrame to long format for the two difference columns
    df_long = df.melt(
        id_vars=['Dataset', 'Model'],
        value_vars=['FedscGen - scGen', 'FedscGen-SMPC - scGen'],
        var_name='Method',
        value_name='Accuracy Difference'
    )
    df_long['Dataset'] = df_long['Dataset'].map(dict(zip(DATASETS, DATASETS_ACRONYM)))
    # Update combination_order with new names
    plt.figure(figsize=(15, 8))
    ax = sns.boxplot(
        x='Dataset',
        y='Accuracy Difference',
        hue='Combination',
        data=df_long,
        palette='Set2',
    )
    set_fontsize(ax, 'Accuracy Difference', 14, 12)
    # Customize the plot
    plt.title('Cell type classification accuracy difference of corrected data (Using FedscGen vs scGen)', fontsize=14)
    plt.xticks(rotation=45, ha='right')
    plt.axhline(0, color='gray', linestyle='--', linewidth=1)  # Add a reference line at 0
    plt.legend(title='', bbox_to_anchor=(1.05, 1), loc='upper left')
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    plt.savefig(f'{plot_dir}/classification_accuracy_difference.png', dpi=300)


def get_classification_stats(data_dir):
    """Reads classification accuracy files, aligns runs, and computes statistics.

    Returns a DataFrame with columns: Dataset, Model, Run, Approach, Accuracy, AUC
    """
    results = []
    n_clients_map = {"HumanPancreas": 5, "CellLine": 3}

    for approach in ["scGen", "FedscGen", "FedscGen-SMPC"]:
        for dataset in DATASETS:
            if dataset == "MouseBrain":
                continue
            files_dir = os.path.join(data_dir, approach.lower(), dataset)

            if approach == "scGen":
                files_path = f"{files_dir}/all/classification/classification_acc_*.csv"
                files = glob.glob(files_path)
            elif approach == "FedscGen":
                n_clients = n_clients_map.get(dataset, 2)
                files_path = f"{files_dir}/all/BO0-C{n_clients}/classification/classification_acc_*.csv"
                files = glob.glob(files_path)
            else:  # FedscGen-SMPC
                files_path = f"{files_dir}/classification/classification_acc_*.csv"
                files = glob.glob(files_path)

            if not files:
                raise FileNotFoundError(f"No classification accuracy files found for in {files_path}")

            for file in files:
                try:
                    fold = int(file.split("_")[-1].split(".csv")[0])  # Extract run number
                    df = pd.read_csv(file)

                    # Get best epoch based on Accuracy
                    best_epoch = df.loc[df["Accuracy"].idxmax()]
                    max_acc = best_epoch["Accuracy"]
                    max_auc = best_epoch["AUC"]

                    results.append({
                        "Dataset": dataset,
                        "Fold": fold,
                        "Approach": approach,
                        "Accuracy": max_acc,
                        "AUC": max_auc
                    })

                except Exception as e:
                    print(f"[ERROR] Failed to process {file}: {e}")

    return pd.DataFrame(results)


def compute_accuracy_differences(df):
    """Aligns classification stats and computes pairwise accuracy differences from DataFrame."""
    aligned_data = []

    grouped = df.groupby(['Dataset', 'Fold'])

    for (dataset, fold), group in grouped:
        acc_map = {row['Approach']: row['Accuracy'] for _, row in group.iterrows()}

        if all(approach in acc_map for approach in ['scGen', 'FedscGen', 'FedscGen-SMPC']):
            aligned_data.append({
                "Dataset": dataset,
                "Model": "MLP",  # Hardcoded since model info is not included anymore
                "Fold": fold,
                "FedscGen - scGen": acc_map['FedscGen'] - acc_map['scGen'],
                "FedscGen-SMPC - scGen": acc_map['FedscGen-SMPC'] - acc_map['scGen']
            })

    return pd.DataFrame(aligned_data)


def classification_error_bar_plot(data_dir):
    df_stats = get_classification_stats(data_dir)
    accuracy_diff_df = compute_accuracy_differences(df_stats)
    accuracy_diff_df.to_csv(os.path.join(data_dir, "accuracy_diff.csv"), index=False)
    plot_accuracy_diff(accuracy_diff_df, data_dir)


def plot_smpc_wilcoxon_heatmap(data_dir, plot_name="smpc-wilcoxon-heatmap.png"):
    df = read_smpc_wilcoxon_benchmarks(data_dir)
    diff_df = validate_symmetry(df, dict(zip(DATASETS, DATASETS_ACRONYM)), data_dir)
    datasets = df.Dataset.unique()
    corrected_p_values, avg_data = get_corrected_p_values(diff_df, datasets, diff_df.Seed.unique())
    avg_df = pd.DataFrame(avg_data, index=diff_df.Dataset.unique(), columns=diff_df.Metric.unique())
    wilcoxon_diff_heatmap(diff_df, avg_df, corrected_p_values, datasets, f"{data_dir}/{plot_name}")
    df.to_csv(os.path.join(data_dir, "smpc_wilcoxon.csv"), index=False)


def wilcoxon_diff_heatmap(df, avg_df, p_values, datasets, file_path):
    metrics = df.columns
    annotations = avg_df.copy().astype(str)

    for i in range(len(datasets)):
        for j in range(len(metrics)):
            annotations.iloc[i, j] = f"{avg_df.iloc[i, j]:.2f}"  # FIXED HERE

    # Create heatmap
    fig, ax = plt.subplots(1,1,figsize=(16, 14), squeeze=True)
    fig.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15, wspace=0.01, hspace=0.01)
    abs_max = 1  # Max range for colormap
    sns.heatmap(avg_df, annot=annotations, fmt="", cmap='RdBu', center=0, ax=ax, cbar=False,
                annot_kws={"size": 32}, square=True, vmin=-abs_max, vmax=abs_max)
    # Set title and labels
    ax.set_title("")
    ax.grid(False)
    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=30, rotation=0, ha='right')
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=30, rotation=45, ha='right')
    # Add colorbar
    norm = colors.Normalize(vmin=-abs_max, vmax=abs_max)
    mappable = plt.cm.ScalarMappable(cmap='RdBu', norm=norm)
    cbar_ax = fig.add_axes([0.81, 0.25, 0.02, 0.5])
    cbar = plt.colorbar(mappable, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=24)
    # Manually add stars above numbers
    for i in range(len(datasets)):
        for j in range(len(metrics)):
            p = p_values[i, j]
            stars = ""
            if p < 0.01:
                stars = "***"
            elif p < 0.05:
                stars = "**"
            elif p < 0.1:
                stars = "*"

            if stars:
                ax.text(j + 0.5, i + 0.25, stars, fontsize=28, ha='center', va='bottom', color='black',
                        fontweight='bold')
    plt.savefig(file_path, dpi=300)
    plt.close()

def get_corrected_p_values(df, datasets, seeds):
    """
    Computes Wilcoxon signed-rank test for each dataset-metric pair across all seeds.
    Applies FDR correction for multiple testing.
    """
    df_wide = df.pivot(index=['Seed', 'Dataset'], columns='Metric', values='Difference')
    df_wide.reset_index(inplace=True)

    # Optional: Sort by Dataset and Seed if needed
    df_wide = df_wide.sort_values(by=['Dataset', 'Seed']).reset_index(drop=True)
    df = df_wide.drop(columns=['Seed', 'Dataset'])
    num_datasets = len(datasets)
    num_metrics = len(df.columns)
    num_seeds = len(seeds)

    # Reshape data to handle seeds explicitly
    avg_data = np.zeros((num_datasets, num_metrics))  # Store mean performance differences
    p_values = np.full((num_datasets, num_metrics), np.nan)  # Initialize with NaN

    for i, dataset in enumerate(datasets):
        start_idx = i * num_seeds
        end_idx = (i + 1) * num_seeds
        dataset_data = df.iloc[start_idx:end_idx]  # Extract data for this dataset

        for j, metric in enumerate(df.columns):
            metric_values = dataset_data[metric].values

            if np.std(metric_values) > 0:  # Ensure variation exists
                # Wilcoxon signed-rank test across all seeds for this dataset-metric pair
                stat, p = wilcoxon(metric_values)
                p_values[i, j] = p  # Store p-value

            # Compute average performance difference
            avg_data[i, j] = np.mean(metric_values)

    # Multiple testing correction (FDR)
    valid_p_values = p_values[~np.isnan(p_values)]  # Filter out NaNs
    corrected_p = multipletests(valid_p_values, method='fdr_bh')[1]

    # Store corrected p-values
    p_values_corrected = np.full_like(p_values, np.nan)
    p_values_corrected[~np.isnan(p_values)] = corrected_p

    return p_values_corrected, avg_data


def read_smpc_wilcoxon_benchmarks(data_dir, filename="benchmark_metrics.csv"):
    """
    Parameters
    ----------
    data_dir: str
        Assuming data_dir has two subfolders: "scgen" and "fedscgen-smpc"
    Returns
    -------
    """
    def read_benchmarks(approach):
        file_path = os.path.join(data_dir, approach, filename)
        if os.path.exists(file_path):
            df = pd.read_csv(file_path)
        else:
            raise FileNotFoundError(file_path)
        df = df[df.Batch.isna() & df.BatchOut.isna() & (df.Inclusion == "all") & (df.Dataset != "MouseBrain")]
        df.drop(columns=["Inclusion", "Epoch", "Round", "Batch", "BatchOut", "N_Clients", "File"], inplace=True)
        return df
    fedscgen_smpc_df = read_benchmarks("fedscgen-smpc")
    scgen_df = read_benchmarks("scgen")
    def sanity_check(df, name):
        print(f"\n--- Sanity check for {name} ---")
        dataset_group_sizes = df.groupby(["Dataset", "Seed"]).size().unstack(fill_value=0)
        valid_seeds = dataset_group_sizes.columns[
            (dataset_group_sizes != 0).all(axis=0)
        ].tolist()
        print(f"Dropped seeds: {set(dataset_group_sizes.columns) - set(valid_seeds)}")
        # Filter dataframe to keep only rows with shared seeds
        df_filtered = df[df.Seed.isin(valid_seeds)].copy()
        return df_filtered
    fedscgen_smpc_df = sanity_check(fedscgen_smpc_df, "fedscgen-smpc")
    scgen_df = sanity_check(scgen_df, "scgen")
    sort_cols = ["Dataset", "Seed"]
    fedscgen_smpc_df = fedscgen_smpc_df.sort_values(by=sort_cols).reset_index(drop=True)
    scgen_df = scgen_df.sort_values(by=sort_cols).reset_index(drop=True)
    assert fedscgen_smpc_df["Dataset"].equals(scgen_df["Dataset"]), "❌ Dataset values do not match!"
    assert fedscgen_smpc_df["Seed"].equals(scgen_df["Seed"]), "❌ Seed values do not match!"
    combined_df = pd.concat([fedscgen_smpc_df, scgen_df], ignore_index=True)
    return combined_df

def validate_symmetry(df, ds_map, plot_dir):
    df_pivot = df.pivot(index=['Seed', 'Dataset'], columns='Approach')
    metrics = [c for c in df.columns if c not in ['Seed', 'Dataset', 'Approach']]
    df_diff = df_pivot.xs('FedscGen-SMPC', axis=1, level=1)[metrics] - df_pivot.xs('scGen', axis=1, level=1)[metrics]
    df_diff.reset_index(inplace=True)
    df = df_diff.melt(id_vars=['Seed', 'Dataset'], var_name='Metric', value_name='Difference')
    df = df[df.Metric != "ILF1"]
    df = df[df.Metric != "GC"]
    if 'Unnamed: 0' in df.columns:
        df = df.drop(columns=['Unnamed: 0'])
    df.Dataset = df.Dataset.apply(lambda x: ds_map[x])
    datasets = df.Dataset.unique()
    metrics =  df.Metric.unique()
    seeds = df.Seed.unique()
    diff_df = df # compute_differences(df, datasets, metrics, seeds)
    plot_all_distributions(diff_df, save_path=f"{plot_dir}/symmetry_grid.png")
    symmetry_df = calculate_symmetry_metrics(diff_df, datasets, metrics)
    symmetry_df.to_csv(f"{plot_dir}/symmetry_metrics.csv", index=False)
    # check_symmetry_holds(symmetry_df, skewness_threshold=0.5, bowley_threshold=0.3)
    check_symmetry_by_skewness(symmetry_df, skewness_threshold=1)
    return diff_df

def plot_all_distributions(diff_df, save_path="diff_grid_metrics_rows_with_inset.png"):
    metrics = sorted(diff_df["Metric"].unique())
    datasets = sorted(diff_df["Dataset"].unique())
    n_rows = len(metrics)
    n_cols = len(datasets)

    # Single subplot per dataset-metric pair (QQ with inset histogram)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3.5 * n_rows))

    # Handle case where there's only one row or column
    if n_rows == 1:
        axes = np.array([axes])
    if n_cols == 1:
        axes = axes.reshape(-1, 1)

    for i, metric in enumerate(metrics):
        for j, dataset in enumerate(datasets):
            diffs = diff_df[
                (diff_df["Dataset"] == dataset) & (diff_df["Metric"] == metric)
                ]["Difference"].values

            # QQ plot
            ax = axes[i][j]
            stats.probplot(diffs, dist="norm", plot=ax)

            # Style the QQ plot
            ax.get_lines()[0].set_marker('o')
            ax.get_lines()[0].set_color('blue')
            ax.get_lines()[1].set_color('red')
            ax.get_lines()[1].set_linewidth(2)

            # Add inset for histogram + KDE (upper left, like your example)
            ax_inset = inset_axes(ax, width="40%", height="40%", loc="upper left")
            sns.histplot(diffs, kde=True, bins=5, ax=ax_inset, color="lightblue", edgecolor="gray")
            ax_inset.set_title("", fontsize=8)
            ax_inset.set_xlabel("", fontsize=6)
            ax_inset.set_ylabel("", fontsize=6)
            ax_inset.tick_params(labelsize=6)
            xlim = ax_inset.get_xlim()
            x_ticks = [xlim[0], xlim[1]]
            ax_inset.set_xticks(x_ticks)
            ax_inset.set_xticklabels([f"{x:.2f}" for x in x_ticks], fontsize=8)

            ylim = ax_inset.get_ylim()
            y_ticks = [1, int(np.floor(max(ylim)))]
            ax_inset.set_yticks(y_ticks)
            ax_inset.set_yticklabels(y_ticks, fontsize=12)
            ax_inset.yaxis.tick_right()
            ax_inset.spines["left"].set_visible(False)
            ax_inset.spines["top"].set_visible(False)
            ax_inset.tick_params(axis="y", which="both", left=False)
            ax_inset.tick_params(axis="x", which="both", top=False)
            ax.set_title("")
            # Labels
            # set off the y-axis labels
            if i == 0:
                ax.set_title(f"{dataset}", fontsize=22)
            if j == 0:
                ax.set_ylabel(f"{metric}", fontsize=22)
                ylim = ax.get_ylim()
                y_ticks = np.linspace(ylim[0], ylim[1], 4)
                ax.set_yticks(y_ticks)
                ax.set_yticklabels([f"{y:.2f}" for y in y_ticks], fontsize=14)
            else:
                ax.set_ylabel("")
                ax.get_yaxis().set_visible(False)

            ax.set_xlabel("")
            if i == n_rows - 1:
                xlim = ax.get_xlim()
                xmin, xmax = xlim
                x_ticks = [round(xmin, 2), 0, round(xmax, 2)]
                ax.set_xticks(x_ticks)
                ax.set_xticklabels([f"{x:.2f}" for x in x_ticks], fontsize=14)
            else:
                ax.get_xaxis().set_visible(False)


    plt.tight_layout(pad=.1, h_pad=0.1, w_pad=0.1)
    plt.savefig(save_path, dpi=300)
    plt.close()


# Function to calculate skewness and Bowley’s skewness
def calculate_symmetry_metrics(diff_df, datasets, metrics):
    symmetry_results = []
    for dataset in datasets:
        for metric in metrics:
            diffs = diff_df[(diff_df["Dataset"] == dataset) & (diff_df["Metric"] == metric)]["Difference"].values
            # skewness = stats.skew(diffs)
            skewness = median_skewness(diffs)
            q1, q2, q3 = np.percentile(diffs, [25, 50, 75])
            bowley = (q3 + q1 - 2 * q2) / (q3 - q1) if (q3 - q1) != 0 else np.nan
            symmetry_results.append({
                "Dataset": dataset,
                "Metric": metric,
                "Skewness": skewness,
                "Bowley_Skewness": bowley
            })
    return pd.DataFrame(symmetry_results)


def check_symmetry_holds(symmetry_df, skewness_threshold, bowley_threshold):
    # Filter tests where both skewness and Bowley's are within bounds
    symmetric_tests = symmetry_df[
        (symmetry_df["Skewness"].abs() <= skewness_threshold) &
        (symmetry_df["Bowley_Skewness"].abs() <= bowley_threshold)
        ]

    # Total tests
    total_tests = len(symmetry_df)
    symmetric_count = len(symmetric_tests)

    # Print summary
    print(f"Total tests: {total_tests}")
    print(
        f"Tests satisfying symmetry assumption (|Skewness| <= {skewness_threshold}, |Bowley| <= {bowley_threshold}): {symmetric_count}")

    # Print symmetric test names in Metric-Dataset format
    print("\nSymmetric Tests (Metric-Dataset):")
    for index, row in symmetric_tests.iterrows():
        test_name = f"{row['Metric']}-{row['Dataset']}"
        print(f"{test_name} (Skewness: {row['Skewness']:.3f}, Bowley: {row['Bowley_Skewness']:.3f})")

    # Print failing tests (if any)
    asymmetric_tests = symmetry_df[
        ~((symmetry_df["Skewness"].abs() <= skewness_threshold) &
          (symmetry_df["Bowley_Skewness"].abs() <= bowley_threshold))
    ]
    if not asymmetric_tests.empty:
        print("\nTests failing symmetry assumption (Metric-Dataset):")
        for index, row in asymmetric_tests.iterrows():
            test_name = f"{row['Metric']}-{row['Dataset']}"
            print(f"{test_name} (Skewness: {row['Skewness']:.3f}, Bowley: {row['Bowley_Skewness']:.3f})")
    symmetric_tests.to_csv("symmetric_tests.csv", index=False)
    asymmetric_tests.to_csv("asymmetric_tests.csv", index=False)

def median_skewness(diffs):
    median = np.median(diffs)
    n = len(diffs)
    mean_dev = diffs - median
    skewness_median = (np.mean(mean_dev**3)) / (np.mean(mean_dev**2)**1.5)
    return skewness_median

def check_symmetry_by_skewness(symmetry_df, skewness_threshold):
    # Filter tests where skewness is within bounds
    symmetric_tests = symmetry_df[
        (symmetry_df["Skewness"].abs() <= skewness_threshold)
    ]

    # Total tests
    total_tests = len(symmetry_df)
    symmetric_count = len(symmetric_tests)

    # Print summary
    print(f"Total tests: {total_tests}")
    print(f"Tests satisfying symmetry assumption (|Skewness| <= {skewness_threshold}): {symmetric_count}")

    # Print symmetric test names in Metric-Dataset format
    print("\nSymmetric Tests (Metric-Dataset):")
    for index, row in symmetric_tests.iterrows():
        test_name = f"{row['Metric']}-{row['Dataset']}"
        print(f"{test_name} (Skewness: {row['Skewness']:.3f}, Bowley: {row['Bowley_Skewness']:.3f})")

    # Print failing tests
    asymmetric_tests = symmetry_df[
        (symmetry_df["Skewness"].abs() > skewness_threshold)
    ]
    if not asymmetric_tests.empty:
        print("\nTests failing symmetry assumption (Metric-Dataset):")
        for index, row in asymmetric_tests.iterrows():
            test_name = f"{row['Metric']}-{row['Dataset']}"
            print(f"{test_name} (Skewness: {row['Skewness']:.3f}, Bowley: {row['Bowley_Skewness']:.3f})")

if __name__ == '__main__':
    parser = argparse.ArgumentParser("Plot lineplot for the metrics")
    parser.add_argument('--scenario', type=str, choices=['datasets', 'tuning', 'batchout', 'kbet-diff', "smpc-wilcoxon",
                                                         "scenarios", "classification", "lisi", "classification_error_bar"],
                        default='datasets')
    parser.add_argument('--data_dir', type=str, help='data directory')
    parser.add_argument('--bo_metrics_file', type=str,
                        help="Path to the file containing the metrics for all batch out of HP dataset.")
    parser.add_argument('--bo_kbet_dir', type=str,
                        help="Path to the main directory containing the batch out kBET results.")

    parser.add_argument("--all_ds_metrics_file", type=str,
                        help="Path to the file containing the metrics for all datasets.")

    parser.add_argument("--hp_kbet_file", type=str, help="Path to the file containing the kBET results for HP dataset.")
    parser.add_argument("--output_dir", type=str, help="Path to the output directory.")
    args = parser.parse_args()
    if args.scenario == "tuning":
        read_tuning_res(args.data_dir, True)
    elif args.scenario == "kbet-diff":
        read_kbet(args.data_dir)
    elif args.scenario == "batchout":
        read_batchout(args.bo_metrics_file, args.hp_kbet_file, args.bo_kbet_dir, args.all_ds_metrics_file,
                      args.output_dir)
    elif args.scenario == "scenarios":
        read_scenarios_metrics(args.data_dir)
    elif args.scenario == "datasets":
        read_datasets_metrics(args.data_dir)
    elif args.scenario == "classification":
        read_classification(args.data_dir)
    elif args.scenario == "classification_error_bar":
        classification_error_bar_plot(args.data_dir)
    elif args.scenario == "lisi":
        read_lisi(args.data_dir)
    elif args.scenario == "smpc-wilcoxon":
        plot_smpc_wilcoxon_heatmap(args.data_dir)
