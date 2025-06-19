import copy
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from functools import partial
import matplotlib.lines as mlines


def generate_palette(unique_celltypes):
    """
    Build a large palette:
     - First 20 colors from 'tab20'
     - If > 20, next up to 20 from 'tab20b'
     - If > 40, sample the remainder from 'gist_ncar'
    Args:
        unique_celltypes: list of category names
    Returns:
        dict mapping each category to an (r, g, b, a) tuple
    """
    n_cats = len(unique_celltypes)
    palette_ = {}

    if n_cats <= 20:
        base = plt.get_cmap("tab20", 20)
        for i, cat in enumerate(unique_celltypes):
            palette_[cat] = base(i)
    elif n_cats <= 40:
        base1 = plt.get_cmap("tab20", 20)
        base2 = plt.get_cmap("tab20b", 20)
        for i, cat in enumerate(unique_celltypes):
            if i < 20:
                palette_[cat] = base1(i)
            else:
                palette_[cat] = base2(i - 20)
    else:
        base1 = plt.get_cmap("tab20", 20)
        base2 = plt.get_cmap("tab20b", 20)
        # For the remainder, sample evenly from gist_ncar
        n_extra = n_cats - 40
        base3 = plt.get_cmap("gist_ncar", n_extra)
        for i, cat in enumerate(unique_celltypes):
            if i < 20:
                palette_[cat] = base1(i)
            elif i < 40:
                palette_[cat] = base2(i - 20)
            else:
                palette_[cat] = base3(i - 40)
    return palette_

def umap_plot(data, wspace, plt_name, batch_key, cell_label_key, use_rep=None):
    """
    Generate UMAP plots for batch and cell type annotations using a custom palette.

    Args:
        data: AnnData object
        wspace: Width space between subplots
        plt_name: Name for saving the plot
        batch_key: Key in data.obs for batch annotations
        cell_label_key: Key in data.obs for cell type annotations
        use_rep: Representation to use for neighbors (optional)

    Returns:
        None (saves UMAP plot)
    """
    # Extract unique categories
    batch_categories = data.obs[batch_key].cat.categories.tolist()
    cell_categories = data.obs[cell_label_key].cat.categories.tolist()

    # Generate color palettes
    batch_color_dict = generate_palette(batch_categories)
    cell_color_dict = generate_palette(cell_categories)

    # Compute neighbors and UMAP
    if use_rep:
        sc.pp.neighbors(data, use_rep=use_rep)
    else:
        sc.pp.neighbors(data)
    sc.tl.umap(data)

    # Combine color dictionaries
    combined_color_dict = {**batch_color_dict, **cell_color_dict}

    # Generate UMAP plot
    sc.pl.umap(
        data,
        color=[batch_key, cell_label_key],
        frameon=False,
        wspace=wspace,
        save=plt_name,
        palette=combined_color_dict
    )

def plot_all_umaps(uncorrected, corrected, batch_key, cell_key, umap_directory):
    plot = partial(umap_plot, batch_key=batch_key, cell_label_key=cell_key)
    set_uamp_conf(umap_directory)
    plot(uncorrected, wspace=0.6, plt_name="dataset.png")
    plot(corrected, wspace=0.6, plt_name="source_corrected.png")


def set_uamp_conf(umap_directory):
    sc.settings.figdir = umap_directory
    sc.settings.set_figure_params(dpi=200, frameon=False)
    sc.set_figure_params(dpi=200)
    sc.set_figure_params(figsize=(4, 4))


def single_plot(adata, batch_key, cell_key, umap_directory, plot_name):
    set_uamp_conf(umap_directory)
    umap_plot(adata, batch_key=batch_key, cell_label_key=cell_key, plt_name=plot_name, wspace=0.6)


def plot_umap_fed_scgen(model, adata, test_adata, test_batches, batches, batch_key, cell_key, output_dir):
    train_batches = copy.deepcopy(batches)
    for batch in test_batches:
        train_batches.remove(batch)
    trns_test_batches = translate(str(test_batches))
    if len(trns_test_batches) == 0:
        trns_test_batches = "BO0"
    # UMAP entire dataset
    corrected_dataset = model.batch_removal(adata, batch_key=batch_key, cell_label_key=cell_key,
                                            return_latent=True)
    plot_all_umaps(adata, corrected_dataset, batch_key, cell_key,
                   umap_directory=f"{output_dir}/{trns_test_batches}/dataset")

    # UMAP entire train adata
    train_adata = adata[adata.obs[batch_key].isin(train_batches)].copy()
    corrected_train_adata = model.batch_removal(train_adata,
                                                batch_key=batch_key,
                                                cell_label_key=cell_key,
                                                return_latent=True)
    plot_all_umaps(train_adata, corrected_train_adata, batch_key, cell_key,
                   umap_directory=f"{output_dir}/{trns_test_batches}/train")

    # UMAP entire test adata

    if len(test_batches) > 0:
        corrected_test_adata = model.batch_removal(test_adata,
                                                   batch_key=batch_key,
                                                   cell_label_key=cell_key,
                                                   return_latent=True)
        plot_all_umaps(test_adata, corrected_test_adata, batch_key, cell_key,
                       umap_directory=f"{output_dir}/{trns_test_batches}/test")


def ss_agx_kbet_plot(silhouette_scores, average_gene_expression, kbet_scores, plt_name):
    # Define the metrics and their corresponding data
    metrics_data = {
        'Silhouette Score': silhouette_scores,
        'Average Gene Expression': average_gene_expression,
        'kBET Score': kbet_scores
    }

    # Create subplots
    fig, axs = plt.subplots(nrows=len(metrics_data), figsize=(7, 15))

    # Iterate over the metrics and their corresponding data
    for i, (metric, data) in enumerate(metrics_data.items()):
        # print(i, metric)
        data = data.transpose()
        if metric == 'Average Gene Expression':
            data_source = pd.DataFrame({metric: data.loc['Uncorrected'].tolist(),
                                        'Method': ['Uncorrected'] * len(data.loc['Uncorrected'].tolist())})
            data_corrected = pd.DataFrame({metric: data.loc['Corrected'].tolist(),
                                           'Method': ['Corrected'] * len(data.loc['Corrected'].tolist())})
            data_combined = pd.concat([data_source, data_corrected])
            sns.violinplot(x='Method', y=metric, data=data_combined, ax=axs[i])
        else:
            sns.barplot(x=data.index, y=metric, data=data, ax=axs[i])

        axs[i].set_title(metric)

    plt.tight_layout()
    plt.savefig(plt_name)


def subplot(results, output_dir):
    # Convert the results to a DataFrame
    data = []
    for fold, rounds in results.items():
        for round, scores in rounds.items():
            data.append({
                'Fold': fold,
                'Round': round,
                'Loss': scores['loss'],
                'Accuracy': scores['acc'],
                'AUC': scores['auc'],
            })
    df = pd.DataFrame(data)
    df.to_csv(f"{output_dir}/fed_metrics.csv")
    # Convert the round to an integer for easier plotting
    df['Round'] = df['Round'].str.slice(start=6).astype(int)
    # Create subplots
    fig, axes = plt.subplots(3, 1, figsize=(10, 15))

    # Define the handles and labels list
    handles, labels = [], []

    # Plot the loss, accuracy, and AUC for each fold
    for fold, group in df.groupby('Fold'):
        sns.lineplot(ax=axes[0], x='Round', y='Loss', data=group, linestyle="dashed", legend=False)
        sns.lineplot(ax=axes[1], x='Round', y='Accuracy', data=group, linestyle="dashed", legend=False)
        sns.lineplot(ax=axes[2], x='Round', y='AUC', data=group, linestyle="dashed", legend=False)
        handles.append(mlines.Line2D([], [], color=plt.gca().lines[-1].get_color(), linestyle="dashed"))
        # handles.append(line1.lines[0])
        # Add the handles and labels

    avg = df.groupby('Round').mean()
    sns.lineplot(ax=axes[0], x='Round', y='Loss', data=avg, linestyle="solid", legend=False)
    sns.lineplot(ax=axes[1], x='Round', y='Accuracy', data=avg, linestyle="solid", legend=False)
    sns.lineplot(ax=axes[2], x='Round', y='AUC', data=avg, linestyle="solid", legend=False)
    handles.append(mlines.Line2D([], [], color=plt.gca().lines[-1].get_color(), linestyle="solid"))
    # chars_to_remove = "[]()''"
    # trans = str.maketrans('', '', chars_to_remove)
    labels = [translate(fold) for fold in df.Fold.unique()]
    labels.append("Avg")

    axes[0].set_title('Loss')
    axes[1].set_title('Accuracy')
    axes[2].set_title('AUC by Fold and Round')

    # Define the legend for the whole figure, outside the plots
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.99), ncol=len(handles) // 2)
    plt.savefig(f"{output_dir}/fed_metrics.png")


def translate(s):
    chars_to_remove = "[]()'',"
    trans = str.maketrans(' ', '-', chars_to_remove)
    return s.strip().translate(trans)
