import os
import scanpy as sc
import matplotlib.pyplot as plt
import anndata
import pandas as pd
import matplotlib.lines as mlines
import argparse
from utils import DATASETS, HEX_COLORS, MODEL_COLORS, MINORITY_CLASSES
import sys
from pathlib import Path

# Get the parent directory of the package
parent_dir = str(Path(__file__).resolve().parent.parent)

# Add the parent directory to sys.path if not already included
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

from fedscgen.utils import remove_cell_types, combine_cell_types

MINORITY_CLASSES_DICT = {k: v.split(",") for k, v in zip(DATASETS, MINORITY_CLASSES)}


def gen_color_dict(batch_key, cell_label_key, adata):
    # Define your color map using the specified colors
    if not pd.api.types.is_categorical_dtype(adata.obs[batch_key]):
        batch_labels = adata.obs[batch_key].astype('category').unique().tolist()
    else:
        batch_labels = adata.obs[batch_key].cat.categories.tolist()
    cell_labels = adata.obs[cell_label_key].cat.categories.tolist()
    # Check if there are enough colors for the labels
    assert len(cell_labels) <= len(
        HEX_COLORS), f"Not enough colors for cell type labels ==> {len(cell_labels)} <= {len(HEX_COLORS)}"
    assert len(batch_labels) <= len(
        HEX_COLORS), f"Not enough colors for batch labels ==> {len(batch_labels)} <= {len(HEX_COLORS)}"

    # Create dictionaries mapping labels to colors
    batch_color_dict = {label: HEX_COLORS[i] for i, label in enumerate(batch_labels)}
    cell_color_dict = {label: HEX_COLORS[i] for i, label in enumerate(cell_labels)}
    return batch_color_dict, cell_color_dict


def plot_umap(adata, batch_key, cell_key, dataset_name, inclusion, model, output_dir, batch_color_dict,
              cell_color_dict, plot_legend=False):
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)
    umap_title = f'{dataset_name}-{inclusion}-{model}'
    border_color = MODEL_COLORS[model]

    # Define the output file paths
    batch_output_path = os.path.join(output_dir, f'{umap_title}_batch.png')
    cell_output_path = os.path.join(output_dir, f'{umap_title}_cell.png')

    adata.obs[batch_key] = adata.obs[batch_key].astype('category')
    adata.obs[cell_key] = adata.obs[cell_key].astype('category')

    ind_umap_plot(adata, batch_color_dict, batch_key, batch_output_path, border_color=border_color, linestyle="solid")
    ind_umap_plot(adata, cell_color_dict, cell_key, cell_output_path, border_color=border_color, linestyle="dotted")

    # Save legends in separate files
    # save_legend(batch_color_dict, os.path.join(output_dir, f'{umap_title}_batch_legend.png'))
    # save_legend(cell_color_dict, os.path.join(output_dir, f'{umap_title}_cell_type_legend.png'))
    if model == "Raw" or plot_legend:
        save_combined_legend(batch_color_dict, cell_color_dict, os.path.join(output_dir, f'{umap_title}_legend.png'),
                             dataset_name)


def ind_umap_plot(adata, color_dict, color_key, output_path, border_color='black', border_width=0, linestyle="solid"):
    ax = sc.pl.umap(adata, color=color_key, palette=color_dict, show=False, save=None, legend_loc=None)
    ax.set_title('')
    ax.set_axis_off()
    plt.savefig(output_path, bbox_inches='tight', pad_inches=0.1, dpi=300)
    plt.clf()


def save_legend(color_dict, output_path):
    # Creating an empty plot to steal legend from
    fig, ax = plt.subplots(figsize=(2, len(color_dict) * 0.2))  # Adjust figsize as needed
    ax.axis('off')

    # Creating line objects for each label, just for creating the legend
    lines = [mlines.Line2D([], [], color=color, marker='o', linestyle='', markersize=10)
             for color in color_dict.values()]
    labels = list(color_dict.keys())

    # Creating legend from the line objects
    legend = ax.legend(lines, labels, frameon=False, loc='center', handlelength=1)

    # Adjust legend location and size
    legend.get_frame().set_linewidth(0.0)
    legend.get_frame().set_boxstyle('Square', pad=0)
    fig.tight_layout(rect=[0, 0, 0.75, 1])

    # Saving legend as a separate file
    fig.savefig(output_path, bbox_inches='tight', pad_inches=0, dpi=300)
    plt.close(fig)


from matplotlib.gridspec import GridSpec


def save_combined_legend(batch_color_dict, cell_color_dict, output_path, dataset_name):
    fig = plt.figure(figsize=(4, 10))

    gs = GridSpec(2, 1, height_ratios=[len(batch_color_dict),
                                       len(cell_color_dict)])  # Here, the second subplot is twice the height of the first

    # Create subplots
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    # Hide the axes
    ax1.axis('off')
    ax2.axis('off')

    # Set the title of the dataset at the top of the figure
    fig.suptitle(dataset_name, fontsize=18, fontweight='bold', va='top', ha='left', x=0.05,
                 y=0.98)  # Align to the top-left

    # Creating line objects for batch legend
    batch_lines = [mlines.Line2D([], [], color=color, marker='o', linestyle='', markersize=14)
                   for color in batch_color_dict.values()]
    batch_labels = list(batch_color_dict.keys())
    batch_legend = ax1.legend(batch_lines, batch_labels, frameon=False, loc='upper left', handlelength=1, fontsize=14)
    ax1.set_title('Batch', loc='left', pad=-10, fontsize=16)  # Align to the left

    # Creating line objects for cell type legend
    cell_lines = [mlines.Line2D([], [], color=color, marker='o', linestyle='', markersize=14)
                  for color in cell_color_dict.values()]
    cell_labels = list(cell_color_dict.keys())
    cell_legend = ax2.legend(cell_lines, cell_labels, frameon=False, loc='upper left', handlelength=1, fontsize=14)
    ax2.set_title('Cell type', loc='left', pad=-10, fontsize=16)  # Align to the left

    # Adjust the spacing between the axes and the layout to ensure everything fits without overlapping
    plt.subplots_adjust(hspace=0, wspace=0)  # Adjust vertical space between the axes

    # Save the figure containing both legends and the dataset title
    fig.savefig(output_path, bbox_inches='tight', pad_inches=0.1,
                dpi=300)  # pad_inches adds whitespace around the figure
    plt.close(fig)


def plot_tuning(data_dir, fed_data_dir, raw_data_dir, batch_key, cell_key, epoch, round, output_dir):
    print(f"Plotting UMAP for tuning round {round} and epoch {epoch} ...")
    for dataset_name in DATASETS:
        files = {
            'Raw': os.path.join(raw_data_dir, f'{dataset_name}.h5ad'),
            'ScGen': os.path.join(data_dir, dataset_name, "all", 'corrected.h5ad'),
            'FedScGen': os.path.join(fed_data_dir, dataset_name, f"E{epoch}", f'corrected_{round}.h5ad')
        }
        plot_dataset_umap(batch_key, cell_key, dataset_name, "all", files, output_dir)
        break


def plot_inclusion_scenarios(data_dir, fed_data_dir, raw_data_dir, batch_key, cell_key, output_dir):
    print(f"Plotting UMAP for inclusion scenarios ...")
    inclusions = ['all', 'dropped', 'combined']
    for inclusion in inclusions:
        for dataset_name in DATASETS:
            if inclusion != "all" and dataset_name in ["CellLine", "HumanDendriticCells"]:
                continue
            print(f"Plotting UMAP for {dataset_name} for inclusion {inclusion}...")
            n_clients = 5 if dataset_name == "HumanPancreas" else 3 if dataset_name == "CellLine" else 2
            files = {
                'ScGen': os.path.join(data_dir, dataset_name, inclusion, 'corrected.h5ad'),
                'FedScGen': os.path.join(fed_data_dir, dataset_name, inclusion, f"BO0-C{n_clients}",
                                         "fed_corrected.h5ad"),
                'Raw': os.path.join(raw_data_dir, f'{dataset_name}.h5ad')
            }
            plot_dataset_umap(batch_key, cell_key, dataset_name, inclusion, files, output_dir)


def plot_dataset_umap(batch_key, cell_key, dataset_name, inclusion, files, output_dir, plot_legend=False):
    for key, file in files.items():
        if os.path.exists(file):
            adata = anndata.read_h5ad(file)
            if key == "Raw":
                if inclusion == "dropped":
                    adata = remove_cell_types(adata, MINORITY_CLASSES_DICT[dataset_name], cell_key)
                elif inclusion == "combined":
                    adata = combine_cell_types(adata, MINORITY_CLASSES_DICT[dataset_name], cell_key)
            batch_color_dict, cell_color_dict = gen_color_dict(batch_key, cell_key, adata)
            print(f"Plotting UMAP for {key} file for {dataset_name} ...")
            plot_umap(adata, batch_key, cell_key, dataset_name, inclusion, key, output_dir,
                      batch_color_dict, cell_color_dict, plot_legend=plot_legend)
        else:
            raise FileNotFoundError(f"{key} file for {dataset_name} not found. Skipping... \n {file}")


def plot_batch_out(data_dir, fed_data_dir, batch_key, cell_key, n_batches, output_dir):
    print(f"Plotting UMAP for batch out scenarios ...")
    for b in range(n_batches):
        files = {
            'ScGen': os.path.join(data_dir, f"BO{b}", "corrected.h5ad"),
            'FedScGen': os.path.join(fed_data_dir, f"{b}", "fed_corrected_with_new_studies.h5ad")
        }
        plot_dataset_umap(batch_key, cell_key, f"HumanPancreas-BO{b}", "all", files, output_dir, plot_legend=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate and Plot NMI for a set of adata files.')

    # Add arguments
    parser.add_argument("--data_dir", type=str, help="Path to the main data directory")
    parser.add_argument("--fed_data_dir", type=str, help="Path to the federated data directory")
    parser.add_argument('--raw_data_dir', type=str, help="Path to the raw data directory", default=None)
    parser.add_argument("--output_dir", type=str, required=True, help="Path to the output directory")
    parser.add_argument('--cell_key', help='Cell key name.', default="cell_type")
    parser.add_argument('--batch_key', help='Batch key name.', default="batch")
    parser.add_argument('--round', type=int, default=-1, help="Round number for tuning")
    parser.add_argument('--epoch', type=int, default=-1, help="Epoch number for tuning")
    parser.add_argument('--n_batches', type=int, default=5, help="Number of batches for batch out scenario")
    parser.add_argument('--scenario', type=str, choices=['datasets', 'tuning', 'batchout'], default='datasets')

    # Parse the arguments
    args = parser.parse_args()
    if args.scenario == "datasets":
        plot_inclusion_scenarios(args.data_dir, args.fed_data_dir, args.raw_data_dir, args.batch_key, args.cell_key,
                                   args.output_dir)
    elif args.scenario == "tuning":
        plot_tuning(args.data_dir, args.fed_data_dir, args.raw_data_dir, args.batch_key, args.cell_key, args.epoch,
                    args.round,
                    args.output_dir)
    elif args.scenario == "batchout":
        plot_batch_out(args.data_dir, args.fed_data_dir, args.batch_key, args.cell_key, args.n_batches, args.output_dir)
