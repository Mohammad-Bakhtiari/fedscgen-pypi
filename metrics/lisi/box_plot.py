import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import os
import sys
from pathlib import Path

# Get the parent directory of the package
parent_dir = str(Path(__file__).resolve().parent.parent)

# Add the parent directory to sys.path if not already included
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

from metrics.utils import DATASETS, DATASETS_ACRONYM

# DATASETS = ["HumanDendriticCells",
#             "MouseCellAtlas",
#             "HumanPancreas",
#             "PBMC",
#             "CellLine",
#             "MouseRetina",
#             "MouseBrain",
#             "MouseHematopoieticStemProgenitorCells"
#             ]
# DATASETS_ACCRONYM = ["HDC", "MCA", "HP", "PBMC", "CL", "MR", "MB", "MHSPC"]

color_palette = sns.color_palette("viridis", 3)
markers = ['o', 's', '^', 'p', 'D', 'v', '*', 'H']

def set_fontsize(ax, y_label, font_size, tick_fontsize):
    ax.set_xlabel('', fontsize=font_size)
    ax.set_ylabel(y_label, fontsize=font_size)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=tick_fontsize)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=tick_fontsize)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process LISI results.")
    parser.add_argument('--data_dir', type=str, help="Path to the main directory containing the LISI results.")
    args = parser.parse_args()

    inclusions = ["all", "dropped", "combined"]
    # Loop through each target folder.
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
            csv_path = os.path.join(args.data_dir, dataset, inclusion, f"BO0-C{n_clients}", 'lisi_results.csv')
            if os.path.exists(csv_path):
                df = pd.read_csv(csv_path)
                df.rename(columns={"file_name": "Approach"}, inplace=True)
                df["Approach"] = df["Approach"].apply(lambda x: "FedscGen" if x == "fedscgen" else "scGen" if x== "scgen" else "Raw")
                df["Dataset"] = DATASETS_ACRONYM[ind]
                aggregated_data.append(df)
            else:
                raise FileNotFoundError(f"{csv_path} is not exist")

        # Concatenate all dataframes to get a single dataframe for the current target.
        aggregated_df = pd.concat(aggregated_data, ignore_index=True)
        # order the dataframe raw, scgen, fedscgen
        aggregated_df['Approach'] = pd.Categorical(aggregated_df['Approach'], categories=["Raw", "scGen", "FedscGen"], ordered=True)
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
        plt.savefig(f'{args.data_dir}/lisi_{inclusion}.png', dpi=1000)
        plt.close()
